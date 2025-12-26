#!/bin/bash
#SBATCH --job-name=Prepare_ReadMap_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_ReadMap_Jobs_GL3_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_ReadMap_Jobs_GL3_%j.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
shopt -s nullglob

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/mnt/scratch2/users/3057556/GLESNI/X204SC24084950-Z01-F003/01.RawData/"
JOB_DIR="/users/3057556/jobs/GL3_readmap_jobs"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"
# Change to the working directory
cd "$WORK_DIR" || exit 1

# Define the directory tags for sample directories
DIR_TAGS=("E" "L" "P")

# Loop through each sample directory
for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    if [[ -d "$dir" ]]; then
      for file in "$dir"/*_unmapped_trimmed_R1.fastq*; do
        if [[ ! -f "$file" ]]; then
          continue
        fi

        filename=$(basename "$file")
        reads_file_1="$filename"
        reads_file_2="${reads_file_1/_R1/_R2}"
        sample_name=$(basename "$dir")
        outdir="${sample_name}_readmapped"
        contigs="${sample_name}_unmapped_trimmed_metaspades/contigs.fasta"
        bowtie_db="${sample_name}_readmapped/${sample_name}_bowtie2_db"
        readmap_output_bam_sorted="${sample_name}_readmapped/${sample_name}_readmapped_sorted.bam"
        readmap_contig_summary="${sample_name}_readmapped/${sample_name}_readmapped_contig_summary.txt"
        readmap_cds_summary="${sample_name}_readmapped/${sample_name}_readmapped_cds_summary.txt"


        # existence checks
        if [[ -s "$dir/$readmap_cds_summary" ]]; then
          echo "Skipping $sample_name - readmap_cds_summary results already exist."
          continue
        fi
        if [[ ! -s "$dir/$contigs" ]]; then
          echo "Skipping $sample_name - contigs not found: $dir/$contigs"
          continue
        fi
        if [[ ! -f "$dir/$reads_file_2" ]]; then
          echo "Skipping $sample_name - paired read file not found: $dir/$reads_file_2"
          continue
        fi

        job_script="$JOB_DIR/${sample_name}_readmap_job.sh"

        echo "Directory: $dir"
        echo "Filename: $filename"
        echo "Reads File 1: $reads_file_1"
        echo "Reads File 2: $reads_file_2"
        echo "Sample Name: $sample_name"
        echo "Output Directory: $outdir"
        echo "Contigs File: $contigs"
        echo "Bowtie2 Database: $bowtie_db"
        echo "Sorted BAM: $readmap_output_bam_sorted"
        echo "Contig Summary: $readmap_contig_summary"
        echo "CDS Summary: $readmap_cds_summary"
        echo "Job Script: $job_script"
        echo "----------------------"

        cat <<EOL > "$job_script"
#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --job-name=ReadMap_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=k2-medpri,k2-bioinf,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# load conda (adjust path if needed)
source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

# Create output directory
mkdir -p "$dir/$outdir"

# Build Bowtie2 database
bowtie2-build "$dir/$contigs" "$dir/$bowtie_db" -p 10

# Perform read mapping directly into a sorted BAM file (use a temp prefix)
bowtie2 -x "$dir/$bowtie_db" -1 "$dir/$reads_file_1" -2 "$dir/$reads_file_2" --very-sensitive-local -p 10 | \
samtools view -bS -@ 10 - | \
samtools sort -@ 10 -T "$dir/tmp_${sample_name}_sort" -o "$dir/$readmap_output_bam_sorted"

# Index the sorted BAM file
samtools index "$dir/$readmap_output_bam_sorted" -@ 10

# Generate read mapping summary
total_reads=\$(samtools view -c "$dir/$readmap_output_bam_sorted")
samtools idxstats "$dir/$readmap_output_bam_sorted" | \
awk -v total="\$total_reads" 'BEGIN{print "Contig_Name\\tContig_Length\\tMapped_Reads\\tTotal_Reads\\tProportion"} {mapped=\$3; len=\$2; if (total>0) prop=mapped/total; else prop=0; print \$1"\\t"len"\\t"mapped"\\t"total"\\t"prop}' > "$dir/$readmap_contig_summary"

# Generate CDS summary in one pipeline (no intermediate file)
bedtools intersect -abam "$dir/$readmap_output_bam_sorted" -b "$dir/${sample_name}_pyrodigal/${sample_name}_pyrodigal.gff" -wa -wb -bed | \
awk '{attr=\$NF; if (match(attr,/ID=([^;]+)/,m)){print m[1]} else if (match(attr,/gene=([^;]+)/,m)){print m[1]} else if (match(attr,/Parent=([^;]+)/,m)){print m[1]} else if (match(attr,/Name=([^;]+)/,m)){print m[1]} else {sub(/;.*/,"",attr); print attr}}' | \
 sort | uniq -c > "$dir/$readmap_cds_summary"

EOL

        echo "Submitting job for $sample_name..."
        sbatch "$job_script"

      done
    fi
  done
done
