#!/bin/bash

#SBATCH --job-name=Prepare_MetaBAT2_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_MetaBAT2_Jobs_Sophie.error
#SBATCH --output=/users/3057556/jobs/Prepare_MetaBAT2_Jobs_Sophie.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/mnt/scratch2/users/3057556/SOPHIE/FastQ/"
JOB_DIR="/users/3057556/jobs/metabat2_jobs_sophie"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

if [[ ! -d "$WORK_DIR" ]]; then
    echo "ERROR: WORK_DIR does not exist: $WORK_DIR"
    exit 1
fi

# Change to the working directory
cd "$WORK_DIR" || exit

# Specify data dir tag - First chars present in 'all' data sample dirs to be processed
DIR_TAG="PN0"

# Loop through each sample directory
for dir in "${DIR_TAG}"*; do
    if [[ -d "$dir" ]]; then
        echo "Processing directory: $dir"

        for file in "$dir"/*unmapped_metaspades/contigs.fasta; do
            if [[ -f "$file" ]]; then
                # Store full contigs file path
                contigs="${file}"

                # Look for the *_readmapped directory
                readmapped_dir=$(find "$dir" -type d -name "*_readmapped" | head -n 1)

                # Find a BAM file inside the readmapped directory
                if [[ -n "$readmapped_dir" ]]; then
                    bam=$(find "$readmapped_dir" -type f -name "*_readmapped_sorted.bam" | head -n 1)
                else
                    bam=""
                fi

                # Get sample name
                sample_name=$(basename "$dir")
             
                # Define the job script name
                job_script="$JOB_DIR/${sample_name}_metabat2_job.sh"

                # Echo statements only
                echo "Filename: $filename"
                echo "Contigs File: $contigs"
                echo "BAM File: $bam"
                echo "Sample Name: $sample_name"
                echo "Job Script: $job_script"
                echo "----------------------"

                cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --job-name=MetaBat2_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment handling
module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/metabat2_2.15

# Create output directory
mkdir -p "$dir/$outdir"

# Change to the working directory
cd "$dir/$outdir" || exit

# Run MetaBAT2
runMetaBat.sh -t 10 ../$contigs ../$bam

EOL

                echo "Submitting job for $sample_name..."
                sbatch "$job_script"
            fi
        done
    fi
done
