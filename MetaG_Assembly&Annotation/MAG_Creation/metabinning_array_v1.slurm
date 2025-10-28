#!/bin/bash

#SBATCH --job-name=Prepare_Binning_Jobs
#SBATCH --error=/users/$$/jobs/Prepare_Binning_Jobs_$$.error
#SBATCH --output=/users/$$/jobs/Prepare_Binning_Jobs_$$.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define the path to the working directory
WORK_DIR="/mnt/$$/data"
JOB_DIR="/users/$$/jobs/binning_jobs_$$"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAG="PN"

# Loop through each sample directory
for dir in ./"$DIR_TAG"*; do
  if [[ -d "$dir" ]]; then
    # Find all metaspades output directories
    for metaspades_dir in "$dir"/*_metaspades; do
      if [[ -d "$metaspades_dir" ]]; then
        contigs_file="$metaspades_dir/contigs.fasta"

        # Check if contigs.fasta exists
        if [[ -f "$contigs_file" ]]; then
          # Extract sample name from the metaspades directory
          sample_base=$(basename "$metaspades_dir" _metaspades)

          # Define output directories for binning tools
          metabat2_dir="$dir/${sample_base}_metabat2"
          concoct_dir="$dir/${sample_base}_concoct"
          maxbin2_dir="$dir/${sample_base}_maxbin2"
          metawrap_dir="$dir/${sample_base}_metawrap_refinement"

          # Check if MetaWRAP refinement has already been completed
          if [[ ! -d "$metawrap_dir/metawrap_50_10_bins" ]]; then
            # Define the job script name
            job_script="$JOB_DIR/${sample_base}_binning_job.sh"

            # Get the original fastq file paths for depth calculation
            sample_name_1="${sample_base}_R1.fastq.gz"
            sample_name_2="${sample_base}_R2.fastq.gz"

            # Write the SLURM job script
            cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=30
#SBATCH --mem=350G
#SBATCH --job-name=Binning_${sample_base}
#SBATCH --error=$JOB_DIR/${sample_base}_binning_error.log
#SBATCH --output=$JOB_DIR/${sample_base}_binning_output.log
#SBATCH --partition=k2-medpri
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load the required modules
module load apps/anaconda3/2024.06/bin

echo "Starting binning pipeline for ${sample_base}"
echo "Contigs file: $contigs_file"

# Create output directories
mkdir -p "$metabat2_dir"
mkdir -p "$concoct_dir"
mkdir -p "$maxbin2_dir"
mkdir -p "$metawrap_dir"

# ====================
# Step 1: Run MetaBAT2
# ====================
echo "Running MetaBAT2..."
source activate /mnt/$$/conda-envs/metabat2  # Adjust path as needed

# First, map reads to contigs to generate depth information
# Index the contigs
bowtie2-build "$contigs_file" "$metabat2_dir/${sample_base}_contigs"

# Map reads
bowtie2 -1 "$dir/$sample_name_1" -2 "$dir/$sample_name_2" \
  -x "$metabat2_dir/${sample_base}_contigs" \
  -S "$metabat2_dir/${sample_base}.sam" \
  -p 30

# Convert SAM to sorted BAM
samtools view -bS "$metabat2_dir/${sample_base}.sam" | samtools sort -o "$metabat2_dir/${sample_base}.bam"
samtools index "$metabat2_dir/${sample_base}.bam"

# Generate depth file
jgi_summarize_bam_contig_depths --outputDepth "$metabat2_dir/${sample_base}_depth.txt" "$metabat2_dir/${sample_base}.bam"

# Run MetaBAT2
metabat2 -i "$contigs_file" \
  -a "$metabat2_dir/${sample_base}_depth.txt" \
  -o "$metabat2_dir/${sample_base}_bin" \
  -t 30 \
  -m 1500

conda deactivate

# ====================
# Step 2: Run CONCOCT
# ====================
echo "Running CONCOCT..."
source activate /mnt/$$/conda-envs/concoct  # Adjust path as needed

# Cut up contigs into smaller parts
cut_up_fasta.py "$contigs_file" -c 10000 -o 0 --merge_last -b "$concoct_dir/${sample_base}_contigs_10K.bed" > "$concoct_dir/${sample_base}_contigs_10K.fa"

# Use the BAM file from MetaBAT2 step
# Generate coverage table
concoct_coverage_table.py "$concoct_dir/${sample_base}_contigs_10K.bed" "$metabat2_dir/${sample_base}.bam" > "$concoct_dir/${sample_base}_coverage_table.tsv"

# Run CONCOCT
concoct --composition_file "$concoct_dir/${sample_base}_contigs_10K.fa" \
  --coverage_file "$concoct_dir/${sample_base}_coverage_table.tsv" \
  -b "$concoct_dir/" \
  -t 30

# Merge subcontig clustering into original contig clustering
merge_cutup_clustering.py "$concoct_dir/clustering_gt1000.csv" > "$concoct_dir/${sample_base}_clustering_merged.csv"

# Extract bins
mkdir -p "$concoct_dir/fasta_bins"
extract_fasta_bins.py "$contigs_file" "$concoct_dir/${sample_base}_clustering_merged.csv" --output_path "$concoct_dir/fasta_bins"

conda deactivate

# ====================
# Step 3: Run MaxBin2
# ====================
echo "Running MaxBin2..."
source activate /mnt/$$/conda-envs/maxbin2  # Adjust path as needed

# MaxBin2 needs an abundance file (can use the depth file from MetaBAT2 or create new one)
# Using MetaBAT2 depth file, extract just the depth columns
cut -f1,4 "$metabat2_dir/${sample_base}_depth.txt" | tail -n +2 > "$maxbin2_dir/${sample_base}_abundance.txt"

# Run MaxBin2
run_MaxBin.pl -contig "$contigs_file" \
  -abund "$maxbin2_dir/${sample_base}_abundance.txt" \
  -out "$maxbin2_dir/${sample_base}_bin" \
  -thread 30

conda deactivate

# ====================
# Step 4: Run MetaWRAP Bin_refinement
# ====================
echo "Running MetaWRAP bin refinement..."
source activate /mnt/$$/conda-envs/metawrap  # Adjust path as needed

# MetaWRAP expects bins to be in separate directories with .fa extension
# Prepare bin directories
mkdir -p "$metawrap_dir/metabat2_bins"
mkdir -p "$metawrap_dir/concoct_bins"
mkdir -p "$metawrap_dir/maxbin2_bins"

# Copy MetaBAT2 bins
cp "$metabat2_dir"/*.fa "$metawrap_dir/metabat2_bins/" 2>/dev/null || echo "No MetaBAT2 bins found"

# Copy CONCOCT bins (and rename to .fa if needed)
for bin in "$concoct_dir/fasta_bins"/*.fa; do
  if [[ -f "\$bin" ]]; then
    cp "\$bin" "$metawrap_dir/concoct_bins/"
  fi
done

# Copy MaxBin2 bins (and rename to .fa if needed)
for bin in "$maxbin2_dir"/*.fasta; do
  if [[ -f "\$bin" ]]; then
    bin_name=\$(basename "\$bin" .fasta)
    cp "\$bin" "$metawrap_dir/maxbin2_bins/\${bin_name}.fa"
  fi
done

# Run MetaWRAP bin_refinement module
metawrap bin_refinement \
  -o "$metawrap_dir" \
  -t 30 \
  -A "$metawrap_dir/metabat2_bins" \
  -B "$metawrap_dir/concoct_bins" \
  -C "$metawrap_dir/maxbin2_bins" \
  -c 50 \
  -x 10

conda deactivate

echo "Binning pipeline completed for ${sample_base}"
EOL

            # Submit the job script
            echo "Submitting binning job for $sample_base..."
            sbatch "$job_script"
          else
            echo "Skipping $sample_base: MetaWRAP refinement already completed in $metawrap_dir."
          fi
        else
          echo "Warning: contigs.fasta not found in $metaspades_dir"
        fi
      fi
    done
  fi
done

echo "All binning jobs submitted!"