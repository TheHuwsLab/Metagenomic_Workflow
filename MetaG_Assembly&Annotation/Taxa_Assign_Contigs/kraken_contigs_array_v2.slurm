#!/bin/bash

#SBATCH --job-name=Prepare_Kraken2_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_Kraken2_Jobs_GL6_j5.error
#SBATCH --output=/users/3057556/jobs/Prepare_Kraken2_Jobs_GL6_j%.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/mnt/scratch2/users/3057556/GLESNI/X204SC24084950-Z01-F006/01.RawData"
JOB_DIR="/users/3057556/jobs/GL6_kraken2_contig_jobs"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAGS=("E" "L" "P")

# Define Kraken2 database
DBNAME=/mnt/scratch2/igfs-databases/Kraken2_DBs/k2_pluspfp_20251015/

for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    if [[ -d "$dir" ]]; then
      for file in "$dir"/*_unmapped_metaspades/contigs.fasta; do
        sample_name=$(basename "$dir")
        outdir="${dir}/${sample_name}_kraken2"
        kraken_output="${outdir}/${sample_name}_kraken2.txt"
        report_mpa="${outdir}/${sample_name}_kraken2_report_mpa.txt"
        report="${outdir}/${sample_name}_kraken2_report.txt"

        # Check if output files exist and are non-empty
        if [[ -s "$kraken_output" && -s "$report_mpa" && -s "$report" ]]; then
          echo "Skipping $sample_name - Kraken2 results already exist."
          continue
        fi

        # Create job script
        job_script="$JOB_DIR/${sample_name}_kraken_job.sh"

        cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH --job-name=Kraken2_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=k2-medpri,k2-lowpri,k2-bioinf,k2-himem
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/kraken2_2.17.1

mkdir -p "$outdir"

# Run Kraken2 only if output files are missing or empty
if [[ ! -s "$kraken_output" || ! -s "$report_mpa" || ! -s "$report" ]]; then
  kraken2 --threads 10 --output "$kraken_output" --report "$report_mpa" \\
    --db "$DBNAME" --use-names --use-mpa-style "$file"

  kraken2 --threads 10 --output "$kraken_output" --report "$report" \\
    --db "$DBNAME" --use-names "$file"
else
  echo "Skipping $sample_name - Kraken2 results already exist."
fi
EOL

      # Submit the job
      echo "Submitting job for $sample_name..."
      sbatch "$job_script"
done
    fi
  done
done