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
#WORK_DIR="/mnt/3057556/data"
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
        sample_name_1="${filename}"
        sample_name_2=$(echo "$sample_name_1" | sed 's/R1/R2/g')
        sample_name=$(basename "$dir")  # Extract sample name

        # Define output filenames (Change R1_001.fastq.gz to match your file naming convention)
        sorted_bam=$(echo "$sample_name_1" | sed 's/_R1_001.fastq.gz/_sorted.bam/g')
        unmapped_R1=$(echo "$sample_name_1" | sed 's/_R1_001.fastq.gz/_unmapped_trimmed_R1.fastq.gz/g')
        unmapped_R2=$(echo "$sample_name_2" | sed 's/_R2_001.fastq.gz/_unmapped_trimmed_R2.fastq.gz/g')

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


# Run bowtie2 and pipe the output directly into samtools without writing intermediate files
bowtie2 --threads 20  \
    -x /mnt/scratch2/igfs-databases/reference_genomes/rumen_contamination_genomes/rumen_contamination_bt2_db \
    -1 "$dir/$sample_name_1" -2 "$dir/$sample_name_2" | \
samtools view -@ 10 -b - | \
samtools sort -@ 10 -n - -o "$dir/$sorted_bam"

# Convert unmapped reads to FASTQ
# Compress the output files using pigz
set -o pipefail

samtools fastq -@ 10 -f 12 "$dir/$sorted_bam" \
  -1 >(pigz -p 10 -c > "$dir/${unmapped_R1}") \
  -2 >(pigz -p 10 -c > "$dir/${unmapped_R2}") \
  -0 /dev/null -s /dev/null -n

wait

EOL

      # Submit the job script
      echo "Submitting job for $sample_name..."
      sbatch "$job_script"

done
    fi
  done
done
