#!/bin/bash

#SBATCH --job-name=Prepare_Binning_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_Binning_Jobs_Testing.error
#SBATCH --output=/users/3057556/jobs/Prepare_Binning_Jobs_Testing.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define the path to the working directory
WORK_DIR="/mnt/scratch2/users/3057556/TODO/testing"
JOB_DIR="/users/3057556/jobs/binning_jobs_testing"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAGS=("E" "L" "P" "A")

# Loop through each sample directory
for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    if [[ -d "$dir" ]]; then

      # Find all metaspades output directories
      for metaspades_dir in "$dir"/*_metaspades; do
        ### Search for bam files
        echo "Searching for BAM files..."
        bam_candidates=($(find "$dir" -type f -name "*_readmapped_sorted.bam" 2>/dev/null))
        if [[ ${#bam_candidates[@]} -gt 0 ]]; then
            contig_readmapping_file="${bam_candidates[0]}"
            echo "Using BAM: $contig_readmapping_file"
        else
            echo "No BAM found — will generate via mapping"
            contig_readmapping_file=""
        fi
        ###
        if [[ -d "$metaspades_dir" ]]; then
          contigs_file="$metaspades_dir/contigs.fasta"


          # Check if contigs.fasta exists
          if [[ -f "$contigs_file" ]]; then
            # Extract sample name from the metaspades directory
            sample_base=$(basename "$metaspades_dir" _unmapped_trimmed_paired_metaspades)

            # Define directories
            MAGSCOT_DIR="/mnt/scratch2/igfs-anaconda/conda-envs/magscot/MAGScoT"
            MAG_dir="$dir/${sample_base}_MAGs"
            metabat2_dir="$MAG_dir/${sample_base}_metabat2"
            concoct_dir="$MAG_dir/${sample_base}_concoct"
            maxbin2_dir="$MAG_dir/${sample_base}_maxbin2"
            magscot_dir="$MAG_dir/${sample_base}_magscot_refinement"
            magscot_refined_bins_dir="$magscot_dir/refined_bins"
            checkm2_dir="$MAG_dir/${sample_base}_checkm2"

            # Check if MAGscot refinement has already been completed
            ##if [[ ! -d "$magscot_refined_bins_dir" ]]; then
            # Define the job script name
            job_script="$JOB_DIR/${sample_base}_binning_job.sh"

            # Get the original fastq file paths for depth calculation
            sample_name_1="${sample_base}_unmapped_R1.fastq.gz"
            sample_name_2="${sample_base}_unmapped_R2.fastq.gz"

            # Printcheck: capture and save current variable values for debugging
            capture_file="$JOB_DIR/${sample_base}_capture_check.txt"
            {
              echo "PRINTCHECK for ${sample_base}"
              echo "TAG=${TAG}"
              echo "dir=${dir}"
              echo "metaspades_dir=${metaspades_dir}"
              echo "contig_readmapping_file=${contig_readmapping_file}"
              echo "contigs_file=${contigs_file}"
              echo "sample_base=${sample_base}"
              echo "MAG_dir=${MAG_dir}"
              echo "metabat2_dir=${metabat2_dir}"
              echo "concoct_dir=${concoct_dir}"
              echo "maxbin2_dir=${maxbin2_dir}"
              echo "comebin_dir=${comebin_dir}"
              echo "metawrap_dir=${metawrap_dir}"
              echo "job_script=${job_script}"
              echo "sample_name_1=${sample_name_1}"
              echo "sample_name_2=${sample_name_2}"
              echo ""
              echo "Resolved readmapped dirs:"
              ls -d "$dir"/*_readmapped 2>/dev/null || echo "  (none)"
              echo ""
              echo "contigs_file exists: $( [[ -f "$contigs_file" ]] && echo yes || echo no )"
              echo "metaspades_dir is dir: $( [[ -d "$metaspades_dir" ]] && echo yes || echo no )"
            } > "$capture_file"
            echo "Saved variable capture to $capture_file"
            cat "$capture_file"

            # Write the SLURM job script
            cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --job-name=Binning_${sample_base}
#SBATCH --error=$JOB_DIR/${sample_base}_binning_error.log
#SBATCH --output=$JOB_DIR/${sample_base}_binning_output.log
#SBATCH --partition=k2-bioinf,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load the required modules
module load apps/anaconda3/2024.06/bin

echo "Starting binning pipeline for ${sample_base}"
echo "Contigs file: $contigs_file"


COMPLETENESS_THRESHOLD=80
CONTAMINATION_THRESHOLD=10

# Create output directories
mkdir -p "$MAG_dir"
mkdir -p "$metabat2_dir"
mkdir -p "$concoct_dir"
mkdir -p "$maxbin2_dir"
mkdir -p "$magscot_refined_bins_dir"
mkdir -p "$checkm2_dir"

if [[ -z "$contig_readmapping_file" || ! -f "$contig_readmapping_file" ]]; then
    echo "ERROR: No valid BAM file available before indexing"
    exit 1
fi


echo "indexing contig read mapping BAM file..."
source activate /mnt/scratch2/igfs-anaconda/conda-envs/samtools_1.21
samtools index "$contig_readmapping_file"


# ====================
# Step 1: Run MetaBAT2
# ====================
echo "Running MetaBAT2..."
source activate /mnt/scratch2/igfs-anaconda/conda-envs/metabat2_2.15 # Adjust path as needed


# Generate depth file
jgi_summarize_bam_contig_depths \
  --outputDepth "$MAG_dir/${sample_base}_depth.txt" \
  "$contig_readmapping_file"

# Run MetaBAT2
metabat2 -i "$contigs_file" \
  -a "$MAG_dir/${sample_base}_depth.txt" \
  -o "$metabat2_dir/${sample_base}_bin" \
  -t 20 \
 -m 1500

source deactivate


# ====================
# Step 2: Run CONCOCT
# ====================
echo "Running CONCOCT..."
source activate /mnt/scratch2/igfs-anaconda/conda-envs/concoct_1.1.0  # Adjust path as needed

# Cut up contigs into smaller parts
cut_up_fasta.py "$contigs_file" -c 10000 -o 0 --merge_last -b "$concoct_dir/${sample_base}_contigs_10K.bed" > "$concoct_dir/${sample_base}_contigs_10K.fa"


# Generate coverage table
concoct_coverage_table.py \
  "$concoct_dir/${sample_base}_contigs_10K.bed" \
  "$contig_readmapping_file" \
  > "$concoct_dir/${sample_base}_coverage_table.tsv"

# Run CONCOCT
concoct --composition_file "$concoct_dir/${sample_base}_contigs_10K.fa" \
  --coverage_file "$concoct_dir/${sample_base}_coverage_table.tsv" \
  -b "$concoct_dir/" \
  -t 20

# Merge subcontig clustering into original contig clustering
merge_cutup_clustering.py "$concoct_dir/clustering_gt1000.csv" > "$concoct_dir/${sample_base}_clustering_merged.csv"

# Extract bins
mkdir -p "$concoct_dir/fasta_bins"
extract_fasta_bins.py "$contigs_file" "$concoct_dir/${sample_base}_clustering_merged.csv" --output_path "$concoct_dir/fasta_bins"

source deactivate

# ====================
# Step 3: Run MaxBin2
# ====================
echo "Running MaxBin2..."
source activate /mnt/scratch2/igfs-anaconda/conda-envs/maxbin2_2.2.7  # Adjust path as needed

# MaxBin2 needs an abundance file (can use the depth file from MetaBAT2 or create new one)
# Using MetaBAT2 depth file, extract just the depth columns
cut -f1,4 "$MAG_dir/${sample_base}_depth.txt" | tail -n +2 > "$maxbin2_dir/${sample_base}_abundance.txt"

# Run MaxBin2
run_MaxBin.pl -contig "$contigs_file" \
  -abund "$maxbin2_dir/${sample_base}_abundance.txt" \
  -out "$maxbin2_dir/${sample_base}_bin" \
  -thread 20

source deactivate

# ==========================
# Step 4: Run MAGScoT
# ==========================

echo "Preparing MAGScoT input mappings..."


# Path to your MAGScoT installation — adjust as needed


source activate /mnt/scratch2/igfs-anaconda/conda-envs/magscot


# -------------------------------------------
# Step 4a: ORF prediction with Prodigal
# -------------------------------------------
echo "Running Prodigal ORF prediction..."

echo "Filtering contigs < 1500bp before Prodigal annotation..."

filtered_contigs="$magscot_dir/${sample_base}_contigs_1500bp.fasta"


awk '
  /^>/ {
    if (seq != "" && length(seq) >= 1500) print header"\n"seq
    header=\$0
    seq=""
  }
  !/^>/ { seq = seq \$0 }
  END   { if (length(seq) >= 1500) print header"\n"seq }
' "$contigs_file" > "\$filtered_contigs"


echo "Contigs before filtering: $(grep -c "^>" "$contigs_file")"
echo "Contigs after filtering:  \$(grep -c '^>' \$filtered_contigs)"

if [[ ! -s "\$filtered_contigs" ]]; then
  echo "ERROR: No contigs passed the 1500bp filter — cannot run Prodigal."
  exit 1
fi

##echo "Running Prodigal ORF prediction on filtered contigs..."
prodigal \
  -i "\$filtered_contigs" \
  -p meta \
  -a "$magscot_dir/${sample_base}_prodigal.faa" \
  -d "$magscot_dir/${sample_base}_prodigal.ffn" \
  -o "$magscot_dir/${sample_base}_prodigal.out"

# -------------------------------------------
# Step 4b: HMMer marker gene annotation
# (GTDBtk r207 TIGRfam + Pfam marker HMMs)
# -------------------------------------------
echo "Running HMMer against GTDBtk r207 marker genes..."
hmmsearch \
  -o "$magscot_dir/${sample_base}_hmm_tigrfam.out" \
  --tblout "$magscot_dir/${sample_base}_hmm_tigrfam_hits.out" \
  --noali --notextw --cut_nc \
  --cpu 20 \
  "$MAGSCOT_DIR/hmm/gtdbtk_rel207_tigrfam.hmm" \
  "$magscot_dir/${sample_base}_prodigal.faa"

hmmsearch \
  -o "$magscot_dir/${sample_base}_hmm_pfam.out" \
  --tblout "$magscot_dir/${sample_base}_hmm_pfam_hits.out" \
  --noali --notextw --cut_nc \
  --cpu 20 \
  "$MAGSCOT_DIR/hmm/gtdbtk_rel207_Pfam-A.hmm" \
  "$magscot_dir/${sample_base}_prodigal.faa"

# Parse HMMer output into 3-column format: gene_id, marker, e-value
grep -v "^#" "$magscot_dir/${sample_base}_hmm_tigrfam_hits.out" \
 | awk '{print \$1"\t"\$3"\t"\$5}' \
  > "$magscot_dir/${sample_base}.tigr"

grep -v "^#" "$magscot_dir/${sample_base}_hmm_pfam_hits.out" \
  | awk '{print \$1"\t"\$4"\t"\$5}' \
  > "$magscot_dir/${sample_base}.pfam"

cat "$magscot_dir/${sample_base}.pfam" \
    "$magscot_dir/${sample_base}.tigr" \
    > "$magscot_dir/${sample_base}.hmm"

# -------------------------------------------
# Step 4c: Build 3-column contig-to-bin TSV
# Format: BIN <tab> CONTIG <tab> BINNER
# -------------------------------------------
echo "Building contig-to-bin mapping TSV..."

# MetaBAT2
for bin in "$metabat2_dir"/*.fa; do
  [[ -f "\$bin" ]] || continue
  bin_name=\$(basename "\$bin" .fa)
  grep "^>" "\$bin" | sed 's/>//' | awk -v b="metabat2_\${bin_name}" '{print b"\t"\$1"\tmetabat2"}'
done > "$magscot_dir/${sample_base}_metabat2.contigs_to_bin.tsv"

# CONCOCT
for bin in "$concoct_dir/fasta_bins"/*.fa; do
  [[ -f "\$bin" ]] || continue
  bin_name=\$(basename "\$bin" .fa)
  grep "^>" "\$bin" | sed 's/>//' | awk -v b="concoct_\${bin_name}" '{print b"\t"\$1"\tconcoct"}'
done > "$magscot_dir/${sample_base}_concoct.contigs_to_bin.tsv"

# MaxBin2
for bin in "$maxbin2_dir"/*.fasta; do
  [[ -f "\$bin" ]] || continue
  bin_name=\$(basename "\$bin" .fasta)
  grep "^>" "\$bin" | sed 's/>//' | awk -v b="maxbin2_\${bin_name}" '{print b"\t"\$1"\tmaxbin2"}'
done > "$magscot_dir/${sample_base}_maxbin2.contigs_to_bin.tsv"

# Concatenate all into final input TSV (no header)
magscot_input_tsv="$magscot_dir/${sample_base}_all.contigs_to_bin.tsv"
cat \
  "$magscot_dir/${sample_base}_metabat2.contigs_to_bin.tsv" \
  "$magscot_dir/${sample_base}_concoct.contigs_to_bin.tsv" \
  "$magscot_dir/${sample_base}_maxbin2.contigs_to_bin.tsv" \
  > "\$magscot_input_tsv"



if [[ ! -s "\$magscot_input_tsv" ]]; then
  echo "ERROR: MAGScoT input TSV is empty — no bins produced by any binner. Exiting."
  exit 1
fi

# -------------------------------------------
# Step 4d: Run MAGScoT
# -------------------------------------------
echo "Running MAGScoT refinement..."

Rscript "$MAGSCOT_DIR/MAGScoT.R" \
  -i "\$magscot_input_tsv" \
  --hmm "$magscot_dir/${sample_base}.hmm" \
  -o "$magscot_dir/${sample_base}_magscot"


source deactivate

# -------------------------------------------
# Step 4e: Extract refined bin FASTAs
# MAGScoT outputs <prefix>.refined.contig_to_bin.tsv
# with 2 columns: BIN <tab> CONTIG
# -------------------------------------------
magscot_contig_map="$magscot_dir/${sample_base}_magscot.refined.contig_to_bin.out"

if [[ ! -f "\$magscot_contig_map" ]]; then
  echo "ERROR: MAGScoT refined contig map not found at \$magscot_contig_map"
  exit 1
fi

echo "Extracting MAGScoT refined bin FASTAs..."

python3 - <<PYEOF
import os
from collections import defaultdict

fasta_path   = "$contigs_file"
print("FASTA_PATH="+fasta_path)
mapping_path = "\$magscot_contig_map"
print("MAPPING_PATH="+mapping_path)
out_dir      = "$magscot_refined_bins_dir"
print("OUT_DIR="+out_dir)

# Load contig sequences
seqs = {}
current = None
with open(fasta_path) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            current = line[1:].split()[0]
            seqs[current] = []
        else:
            seqs[current].append(line)

# Load bin assignments — MAGScoT output: BIN <tab> CONTIG
bins = defaultdict(list)
with open(mapping_path) as fh:
    for line in fh:
        parts = line.rstrip().split("\t")
        if len(parts) < 2:
            continue
        bin_id = os.path.basename(parts[0])
        contig_id = parts[1]
        bins[bin_id].append(contig_id)

written = 0
for bin_id, contigs in bins.items():
    out_path = os.path.join(out_dir, f"{bin_id}.fa")
    with open(out_path, "w") as out:
        for contig in contigs:
            if contig in seqs:
                out.write(f">{contig}\n" + "\n".join(seqs[contig]) + "\n")
    written += 1

print(f"Wrote {written} refined bin FASTAs to {out_dir}")
PYEOF

if [[ -z "\$(ls -A "$magscot_refined_bins_dir"/*.fa 2>/dev/null)" ]]; then
  echo "ERROR: No refined bin FASTAs found in $magscot_refined_bins_dir after extraction."
  exit 1
fi

echo "MAGScoT complete. \$(ls "$magscot_refined_bins_dir"/*.fa | wc -l) refined bins produced."

# ============================
# Step 5: Run CheckM2
# ============================
echo "Running CheckM2 on MetaWRAP refined bins..."

checkm2_report="$checkm2_dir/quality_report.tsv"


# Check that MetaWRAP produced bins before proceeding
if [[ -z "\$(ls -A "$magscot_refined_bins_dir"/*.fa 2>/dev/null)" ]]; then
  echo "ERROR: No refined bins found in $refined_bins_dir — skipping CheckM2."
  exit 1
fi

source activate /mnt/scratch2/igfs-anaconda/conda-envs/checkm2_1.1.0

checkm2 predict \
  --threads 20 \
  --input "$magscot_refined_bins_dir" \
  --output-directory "$checkm2_dir" \
  --database_path /mnt/scratch2/igfs-anaconda/conda-dbs/checkm2_1.1.0/CheckM2_database/uniref100.KO.1.dmnd \
  -x fa --force

source deactivate

# Verify the report was produced
if [[ ! -f "\$checkm2_report" ]]; then
  echo "ERROR: CheckM2 quality report not found at \$checkm2_report — cannot filter MAGs."
  exit 1
fi

echo "CheckM2 complete. Filtering MAGs (completeness >= \${COMPLETENESS_THRESHOLD}%, contamination <= \${CONTAMINATION_THRESHOLD}%)..."

# ============================
# Step 6: Filter MAGs by quality
# ============================
passed=0
failed=0

mkdir -p "$MAG_dir/filtered_MAGs"

while IFS=\$'\t' read -r name completeness contamination rest; do
  if awk -v c="\$completeness" -v x="\$contamination" \
       -v ct="\$COMPLETENESS_THRESHOLD" -v xt="\$CONTAMINATION_THRESHOLD" \
       'BEGIN { exit !(c >= ct && x <= xt) }'; then
    src="$magscot_refined_bins_dir/\${name}.fa"
    if [[ -f "\$src" ]]; then
      cp "\$src" "$MAG_dir/filtered_MAGs/${sample_base}_\${name}.fa"
      echo "  PASSED: \${name} (completeness=\${completeness}%, contamination=\${contamination}%)"
      (( passed++ ))
    else
      echo "  WARNING: Bin \${name} passed thresholds but FASTA not found at \$src"
    fi
  else
    echo "  FAILED: \${name} (completeness=\${completeness}%, contamination=\${contamination}%)"
    (( failed++ ))
  fi
done < <(tail -n +2 "\$checkm2_report")

echo ""
echo "=========================================="
echo "CheckM2 filtering summary for ${sample_base}"
echo "  Completeness threshold : >= \${COMPLETENESS_THRESHOLD}%"
echo "  Contamination threshold: <= \${CONTAMINATION_THRESHOLD}%"
echo "  MAGs passed : \$passed"
echo "  MAGs failed : \$failed"
echo "  Filtered MAGs written to: $MAG_dir/filtered_MAGs"
echo "=========================================="

if [[ \$passed -eq 0 ]]; then
  echo "WARNING: No MAGs passed quality thresholds for ${sample_base}."
fi

echo "Full pipeline completed for ${sample_base}"


EOL

              echo "Submitting binning job for $sample_base..."
              sbatch "$job_script"
##            else
##              echo "Skipping $sample_base: MetaWRAP refinement already completed in $metawrap_dir."
##            fi
          else
            echo "Warning: contigs.fasta not found in $metaspades_dir"
          fi
        fi
      done
    fi
  done
done


echo "All binning jobs submitted!"