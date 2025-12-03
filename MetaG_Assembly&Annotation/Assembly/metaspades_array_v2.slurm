#!/bin/bash

#SBATCH --job-name=Prepare_MetaSpades_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_MetaSpades_Jobs_GL3_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_MetaSpades_Jobs_GL3_%j.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define the path to the working directory
WORK_DIR="/mnt/scratch2/users/3057556/data/"
JOB_DIR="/users/3057556/jobs/metaspades_jobs_GL3"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAGS=("E" "L")

for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    if [[ -d "$dir" ]]; then
      for file in "$dir"/*_unmapped_trimmed_R1.fastq.gz; do
        # Get the base filename of the R1 file
        filename=$(basename "$file")
        sample_name_1="${filename}"
        sample_name_2=$(echo "$sample_name_1" | sed 's/R1/R2/g')
        outdir=$(echo "$sample_name_1" | sed 's/_R1.fastq.gz/_metaspades/g') # Here is the output directory name

        # Check if the 'contigs.fasta' file already exists
        # Full path to the 'contigs.fasta' file
        contigs_file="$dir/$outdir/contigs.fasta"
        if [[ ! -s "$contigs_file" ]]; then # Does not exist or is empty
          # Define the job script name
          job_script="$JOB_DIR/${filename%_unmapped_R1.fastq.gz}_metaspades_job.sh"

          # Write the SLURM job script
          cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=300G
#SBATCH --job-name=MetaSpades_${filename%_unmapped_R1.fastq.gz}
#SBATCH --error=$JOB_DIR/${filename%_unmapped_R1.fastq.gz}_error.log
#SBATCH --output=$JOB_DIR/${filename%_unmapped_R1.fastq.gz}_output.log
#SBATCH --partition=k2-medpri,k2-bioinf,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load the required modules
module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/spades_4.2.0

# Run MetaSPAdes
metaspades.py -1 "$dir/$sample_name_1" -2 "$dir/$sample_name_2" -o "$dir/$outdir" -t 20 -m 300
EOL

        # Submit the job script
        echo "Submitting job for $sample_name_1..."
        sbatch "$job_script"
      else
        echo "Skipping $dir: 'contigs.fasta' already exists in $outdir."
fi
      done
    fi
  done
done
