#!/bin/bash

#SBATCH --job-name=Prepare_GL3_Contam_Jobs
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --error=/users/3057556/jobs/Prepare_Contam_Jobs_GL3_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_Contam_Jobs_GL3_%j.log
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths

WORK_DIR="/mnt/scratch2/users/3057556/data/"
JOB_DIR="/users/3057556/jobs/GL3_contamination_jobs"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAGS=("E" "L" "P")

for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    if [[ -d "$dir" ]]; then
      for file in "$dir"/*R1_001_trimmed_paired.fastq.gz; do # Change this pattern if your files have a different naming convention
        filename=$(basename "$file")
        echo "$filename"
        sample_name_1="$filename"
        # Derive mate filename by replacing the first occurrence of _R1 with _R2
        sample_name_2="${sample_name_1/_R1/_R2}"
        sample_name=$(basename "$dir")  # Extract sample name

        # Derive a stable base by removing _R1/_R2 and everything after (handles various suffixes)
        base="${filename%%_R[12]*}"

        # Define output filenames using the base to avoid overwriting original FASTQ files
        sorted_bam="${base}_sorted.bam"
        unmapped_R1="${base}_unmapped_trimmed_paired_R1.fastq.gz"
        unmapped_R2="${base}_unmapped_trimmed_paired_R2.fastq.gz"

        # Define the job script path
        job_script="$JOB_DIR/${sample_name}_contamination_check_job.sh"

        # Echo statements for debugging/logging
        echo "Directory: $dir"
        echo "Sorted Bamfile: $sorted_bam"
        echo "Job Script: $job_script"
        echo "----------------------"

        # Write the SLURM job script
        cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=ContCheck_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_%j.error
#SBATCH --output=$JOB_DIR/${sample_name}_%j.log
#SBATCH --partition=k2-hipri,k2-medpri,k2-lowpri,k2-bioinf
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules and activate conda environment
module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1


SUCCESS_FLAG="$dir/${sample_name}_decontamination_done.flag"
# ensure parent dir exists and make the trap robust (avoid failures inside trap)
mkdir -p "$(dirname "$SUCCESS_FLAG")" || true
trap 'rc=$?; if [ "${rc:-1}" -eq 0 ]; then mkdir -p "$(dirname "$SUCCESS_FLAG")" >/dev/null 2>&1 || true;
printf "%s\n" "$(date --iso-8601=seconds)" > "$SUCCESS_FLAG" 2>/dev/null || true;
 else rm -f "$SUCCESS_FLAG" 2>/dev/null || true; fi' EXIT

# Verify that fastp completed successfully before continuing
FASTP_FLAG="$dir/fastp.succeed"
if [[ -f "$/FASTP_FLAG" ]]; then
  printf "%s: %s found — proceeding with contamination check for sample %s\n" "$(date --iso-8601=seconds)" "$FASTP_FLAG" "$sample_name"
else
  printf "%s: %s NOT found — aborting contamination check for sample %s\n" "$(date --iso-8601=seconds)" "$FASTP_FLAG" "$sample_name" >&2
  exit 1
fi

# Run bowtie2 and pipe the output directly into samtools without writing intermediate files
bowtie2 --threads 20  \
    -x /mnt/scratch2/igfs-databases/reference_genomes/rumen_contamination_genomes/rumen_contamination_bt2_db \
    -1 "$dir/$sample_name_1" -2 "$dir/$sample_name_2" | \
samtools view -@ 10 -b - | \
samtools sort -@ 10 - -o "$dir/$sorted_bam"

# Convert unmapped reads to FASTQ
# Compress the output files using pigz
#set -o pipefail

#samtools fastq -@ 10 \
#    -f 12 \
#    -F 256 -F 2048 \
#    -1 >(pigz -p 10 -c > "$dir/${unmapped_R1}") \
#    -2 >(pigz -p 10 -c > "$dir/${unmapped_R2}") \
#    -0 /dev/null -s /dev/null -n

# Extract unmapped reads uncompressed first
samtools fastq -@ 10 \
    -f 12 \
    -F 256 -F 2048 \
    -1 "$dir/${unmapped_R1%.gz}" \
    -2 "$dir/${unmapped_R2%.gz}" \
    -0 /dev/null -s /dev/null -n

# Then compress once samtools has fully finished
pigz -p 10 -c "$dir/${unmapped_R1%.gz}" > "$dir/${unmapped_R1}" && rm "$dir/${unmapped_R1%.gz}"
pigz -p 10 -c "$dir/${unmapped_R2%.gz}" > "$dir/${unmapped_R2}" && rm "$dir/${unmapped_R2%.gz}"



EOL

      # Submit the job script
      echo "Submitting job for $sample_name..."
      sbatch "$job_script"

done
    fi
  done
done
