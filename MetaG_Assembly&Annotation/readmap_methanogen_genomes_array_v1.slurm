#!/bin/bash

#SBATCH --job-name=Prepare_ReadMapMethanogenGenomes_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_ReadMapMethanogenGenomes_Jobs_Sophie.error
#SBATCH --output=/users/3057556/jobs/Prepare_ReadMapMethanogenGenomes_Jobs_Sophie.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/mnt/scratch2/users/3057556/SOPHIE/FastQ/"
JOB_DIR="/users/3057556/jobs/readmap_target_genomes_jobs_sophie"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit


# Output dir for all output files 
OUTDIR="27_Methanogens_Readmapped"
mkdir -p $OUTDIR

# Bowtie2 database (genomes we are mapping against)
BT_DB="/mnt/scratch2/users/3057556/databases/methanogen_genomes/27_methanogens/27_methanogens"

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

                readmap_output_bam_sorted="$OUTDIR/${sample_name}_readmapped_sorted.bam"


                # Define the job script name
                job_script="$JOB_DIR/${sample_name}_readmap_job.sh"

                # Echo statements only
                echo "Directory: $dir"
                echo "Filename: $filename"
                echo "Reads File 1: $reads_file_1"
                echo "Reads File 2: $reads_file_2"
                echo "Sample Name: $sample_name"
                echo "Output Directory: $outdir"
                echo "Bowtie2 Database: $BT_DB"
                echo "Sorted BAM: $readmap_output_bam_sorted"
                echo "Job Script: $job_script"
                echo "----------------------"

                # Write the SLURM job script
                cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=ReadMapMethanogenGenomes_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

# Create output directory
 



bowtie2 -x "$BT_DB" -1 "$dir/$reads_file_1" -2 "$dir/$reads_file_2" -p 20 --very-sensitive-local | \
samtools view -bS -@ 20 | \
samtools sort -@ 20 -o "$readmap_output_bam_sorted"

# Index the sorted BAM file
samtools index "$readmap_output_bam_sorted" -@ 20





EOL

                # Submit the job script
                echo "Submitting job for $sample_name..."
                sbatch "$job_script"
            fi
        done
    fi
done
