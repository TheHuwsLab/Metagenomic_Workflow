#!/bin/bash

#SBATCH --job-name=Prepare_HUMAnN_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_HUMAnN_Jobs_$$.error
#SBATCH --output=/users/3057556/jobs/Prepare_HUMAnN_Jobs_$$.txt
#SBATCH --partition=$$
#SBATCH --time=01:00:00
#SBATCH --mail-user=$$
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/$$/"
JOB_DIR="/$$/"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

if [[ ! -d "$WORK_DIR" ]]; then
    echo "ERROR: WORK_DIR does not exist: $WORK_DIR"
    exit 1
fi

# Change to the working directory
cd "$WORK_DIR" || exit

# Specify data dir tag - First chars present in 'all' data sample dirs to be processed
DIR_TAG="PN"

# Loop through each sample directory
for dir in "${DIR_TAG}"*; do
    if [[ -d "$dir" ]]; then
        echo "Processing directory: $dir"
        
        for file in "$dir"/*_unmapped_R1.fastq*; do
            if [[ -f "$file" ]]; then
    # Get sample-specific variables
                filename=$(basename "$file")
                R1_fastq="${filename}"
                R2_fastq=$(echo "$R1_fastq" | sed 's/R1/R2/g')
                R1_fasta=$(echo "$R1_fastq" | sed 's/fastq/fasta/g')
                R2_fasta=$(echo "$R2_fastq" | sed 's/fastq/fasta/g')
                combined_fasta=$(echo "$R1_fasta" | sed 's/R1/Combined/g')
                sample_name=$(basename "$dir")
                outdir="${sample_name}_humann"

                # Define the job script name
                job_script="$JOB_DIR/${sample_name}_humann_job.sh"

                # Echo statements only
                echo "Filename: $filename"
                echo "FASTQ File 1: $R1_fastq"
                echo "FASTQ File 2: $R2_fastq"
                echo "FASTA File 1: $R1_fasta"
                echo "FASTA File 2: $R2_fasta"
                echo "Combined FASTA: $combined_fasta"
                echo "Sample Name: $sample_name"
                echo "Output Directory: $outdir"
                echo "Job Script: $job_script"
                echo "----------------------"

                cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=HUMAnN_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=$$
#SBATCH --nodes=1
#SBATCH --time=100:00:00
#SBATCH --mail-user=$$
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment handling
module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/humann_3.9

# Create output directory
mkdir -p "$dir/$outdir"

# Convert FASTQ to FASTA and combine
zcat "$dir/$R1_fastq" | seqtk seq -A - | pigz -p 20 > "$dir/$R1_fasta"
zcat "$dir/$R2_fastq" | seqtk seq -A - | pigz -p 20 > "$dir/$R2_fasta"
zcat "$dir/$R1_fasta" "$dir/$R2_fasta" | pigz -p 20 > "$dir/$combined_fasta"

# Run HUMAnN
humann --input "$dir/$combined_fasta" --output "$dir/$outdir" --threads 20
EOL

                echo "Submitting job for $sample_name..."
                sbatch "$job_script"
            fi
        done
    fi
done
