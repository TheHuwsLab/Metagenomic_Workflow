#!/bin/bash

#SBATCH --job-name=Prepare_AssemblyReadMap_Jobs
#SBATCH --error=/users/$$/jobs/Prepare_AssemblyReadMap_Jobs_$$.error
#SBATCH --output=/users/$$/jobs/Prepare_AssemblyReadMap_Jobs_$$.txt
#SBATCH --partition=$$
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=$$
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/$$/"
JOB_DIR="/$$/"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

# Loop through each sample directory
for dir in ./*; do
    if [[ -d "$dir" ]]; then
        for file in "$dir"/*_unmapped_R1.fastq*; do
            if [[ -f "$file" ]]; then
                # Get sample-specific variables
                filename=$(basename "$file")
                reads_file_1="${filename}"
                reads_file_2=$(echo "$reads_file_1" | sed 's/R1/R2/g')
                sample_name=$(echo "$file" | sed 's/\.\///;s/\/.*//')
                outdir="${sample_name}_readmapped"
                contigs="${sample_name}_unmapped_metaspades/contigs.fasta"
                bowtie_db="${sample_name}_readmapped/${sample_name}_bowtie2_db"
                readmap_output_sam="${sample_name}_readmapped/${sample_name}_readmapped.sam"
                readmap_output_bam="${sample_name}_readmapped/${sample_name}_readmapped.bam"
                readmap_output_bam_sorted="${sample_name}_readmapped/${sample_name}_readmapped_sorted.bam"
                readmap_contig_summary="${sample_name}_readmapped/${sample_name}_readmapped_contig_summary.txt"
                readmap_cds_summary="${sample_name}_readmapped/${sample_name}_readmapped_cds_summary.txt"

                # Define the job script name
                job_script="$JOB_DIR/${sample_name}_readmap_job.sh"

                # Echo statements only
                echo "Directory: $dir"
                echo "Filename: $filename"
                echo "Reads File 1: $reads_file_1"
                echo "Reads File 2: $reads_file_2"
                echo "Sample Name: $sample_name"
                echo "Output Directory: $outdir"
                echo "Contigs File: $contigs"
                echo "Bowtie2 Database: $bowtie_db"
                echo "Readmap Output (SAM): $readmap_output_sam"
                echo "Readmap Output (BAM): $readmap_output_bam"
                echo "Sorted BAM: $readmap_output_bam_sorted"
                echo "Contig Summary: $readmap_contig_summary"
                echo "CDS Summary: $readmap_cds_summary"
                echo "Job Script: $job_script"
                echo "----------------------"

                # Write the SLURM job script
                cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=30
#SBATCH --mem=80G
#SBATCH --job-name=AssemblyReadMap_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=$$
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=$$
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

# Create output directory
mkdir -p "$dir/$outdir"

# Build Bowtie2 database
bowtie2-build "$dir/$contigs" "$dir/$bowtie_db" -p 30

# Perform read mapping
bowtie2 -x "$dir/$bowtie_db" -1 "$dir/$reads_file_1" -2 "$dir/$reads_file_2" -S "$dir/$readmap_output_sam" -p 30 --very-sensitive-local

# Convert SAM to BAM and sort
samtools view -bS "$dir/$readmap_output_sam" -@ 30 > "$dir/$readmap_output_bam"
samtools sort -o "$dir/$readmap_output_bam_sorted" "$dir/$readmap_output_bam" -@ 30
samtools index "$dir/$readmap_output_bam_sorted" -@ 30

# Generate read mapping summary
total_reads=\$(samtools view -c "$dir/$readmap_output_bam_sorted")
samtools idxstats "$dir/$readmap_output_bam_sorted" | \
awk -v total="\$total_reads" 'BEGIN{print "Contig_Name\tMapped_Reads\tTotal_Reads\tProportion"}{print \$1"\t"\$3"\t"total"\t"\$3/total}' > "$dir/$readmap_contig_summary"

# Generate CDS summary
bedtools intersect -a "$dir/$readmap_output_bam_sorted" -b "${sample_name}/${sample_name}_pyrodigal/${sample_name}_pyrodigal.gff" -wa -wb -bed > "$dir/reads_per_gene_new.txt"
awk '{sub(/;.*/, "", \$NF); print \$NF}' "$dir/reads_per_gene.txt" | sort | uniq -c > "$dir/$readmap_cds_summary"
EOL

                # Submit the job script
                echo "Submitting job for $sample_name..."
                sbatch "$job_script"
            fi
        done
    fi
done

