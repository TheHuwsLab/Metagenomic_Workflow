#!/bin/bash

#SBATCH --job-name=Prepare_Yasmin_Contam_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_Yasmin_Contam_Jobs.error
#SBATCH --output=/users/3057556/jobs/Prepare_Yasmin_Contam_Jobs.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/mnt/scratch2/users/3057556/Yasmin/testing/"
JOB_DIR="/users/3057556/jobs/yasmin_contamination_jobs"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

DIR_TAG="PN"

for dir in ./"$DIR_TAG"*; do
  if [[ -d "$dir" ]]; then
    for file in "$dir"/*R1_001.fastq.gz; do
      filename=$(basename "$file")
      sample_name_1="${filename}"
      sample_name_2=$(echo "$sample_name_1" | sed 's/R1/R2/g')
      sample_name=$(basename "$dir")  # Extract sample name

      # Define output filenames
      unmapped_sorted_bam=$(echo "$sample_name_1" | sed 's/_R1_001.fastq.gz/_sorted_unmapped.bam/g')
      unmapped_R1=$(echo "$sample_name_1" | sed 's/_R1_001.fastq.gz/_mapped_R1.fastq/g')
      unmapped_R2=$(echo "$sample_name_1" | sed 's/_R1_001.fastq.gz/_mapped_R2.fastq/g')

      # Define the job script path
      job_script="$JOB_DIR/${sample_name}_contamination_check_job.sh"

      # Echo statements for debugging/logging
      echo "Directory: $dir"
      echo "Unmapped Sorted Bamfile: $unmapped_sorted_bam" 
      echo "Job Script: $job_script"
      echo "----------------------"

      # Write the SLURM job script
      cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=30
#SBATCH --mem=120G
#SBATCH --job-name=ContCheck_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1



# Run bowtie2 and pipe the output directly into samtools without writing intermediate files
bowtie2 --threads 30 --very-sensitive-local \
    -x /mnt/scratch2/igfs-databases/reference_genomes/combined_cow_sheep_human/Human_Sheep_Cow_Combined_Bowtie2/Human_Sheep_Cow_Combined \
    -1 "$dir/$sample_name_1" -2 "$dir/$sample_name_2" | \
samtools view -b -f 12 - | \
samtools sort -n - -o "$dir/$unmapped_sorted_bam"

# Convert unmapped reads to FASTQ
bedtools bamtofastq -i "$dir/$unmapped_sorted_bam" -fq "$dir/$unmapped_R1" -fq2 "$dir/$unmapped_R2"

EOL

      # Submit the job script
      echo "Submitting job for $sample_name..."
      sbatch "$job_script"

    done
  fi
done
