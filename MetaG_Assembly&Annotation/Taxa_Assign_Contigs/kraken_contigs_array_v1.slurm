#!/bin/bash

#SBATCH --job-name=Prepare_Kraken2_Jobs
#SBATCH --error=/users/$$/jobs/Prepare_Kraken2_Jobs_$$.error
#SBATCH --output=/users/$$/jobs/Prepare_Kraken2_Jobs_$$.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/mnt/$$/data"
JOB_DIR="/users/%%/jobs/%%_kraken2_contig_jobs"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAG="PN"

# Define Kraken2 database
DBNAME=/$$/conda-dbs/kraken2/k2_pluspfp_20240904

for dir in "${DIR_TAG}"*; do
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
#SBATCH --mem=200G
#SBATCH --job-name=Kraken2_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=k2-medpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/$$/conda-envs/kraken2_2.1.3

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
