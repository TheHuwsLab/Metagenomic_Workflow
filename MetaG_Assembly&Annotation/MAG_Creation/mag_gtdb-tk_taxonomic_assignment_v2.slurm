#!/bin/bash

#SBATCH --job-name=Prepare_GTDBTK_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_GTDBTK_Jobs.error
#SBATCH --output=/users/3057556/jobs/Prepare_GTDBTK_Jobs.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=00:10:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define the path to the working directory
WORK_DIR="/mnt/scratch2/users/3057556/TODO/testing"
JOB_DIR="/users/3057556/jobs/gtdbtk_jobs"

# GTDB-TK reference data path (adjust if your site stores it elsewhere)
GTDBTK_DATA_PATH="/mnt/scratch2/igfs-anaconda/conda-dbs/gtdbtk_2.7.1/release232"

# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAGS=("E" "L" "P" "A")

# Loop through each sample directory
for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    if [[ -d "$dir" ]]; then

      for metaspades_dir in "$dir"/*_metaspades; do
        if [[ -d "$metaspades_dir" ]]; then

          # Extract sample name (must match naming convention from metabinning_array_v2.slurm)
          sample_base=$(basename "$metaspades_dir" _unmapped_trimmed_paired_metaspades)

          # Define paths consistent with metabinning_array_v2.slurm
          MAG_dir="$dir/${sample_base}_MAGs"
          filtered_MAGs_dir="$MAG_dir/filtered_MAGs"
          gtdbtk_output_dir="$MAG_dir/${sample_base}_filtered_MAGs_gtdbtk"

          # Only proceed if filtered MAGs exist for this sample
          if [[ ! -d "$filtered_MAGs_dir" ]] || [[ -z "$(ls -A "$filtered_MAGs_dir"/*.fa 2>/dev/null)" ]]; then
            echo "Skipping $sample_base: no filtered MAGs found at $filtered_MAGs_dir"
            continue
          fi

          # Skip if GTDB-TK has already been run for this sample
          if [[ -d "$gtdbtk_output_dir" ]] && \
             [[ -n "$(ls -A "$gtdbtk_output_dir"/classify/*.summary.tsv 2>/dev/null)" ]]; then
            echo "Skipping $sample_base: GTDB-TK classification already completed in $gtdbtk_output_dir"
            continue
          fi

          # Write the individual GTDB-TK SLURM job script
          job_script="$JOB_DIR/${sample_base}_gtdbtk_job.sh"

          cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --job-name=GTDBTK_${sample_base}
#SBATCH --error=$JOB_DIR/${sample_base}_gtdbtk_error.log
#SBATCH --output=$JOB_DIR/${sample_base}_gtdbtk_output.log
#SBATCH --partition=k2-hipri,k2-medpri,k2-bioinf,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load the required modules
module load apps/anaconda3/2024.06/bin

echo "Starting GTDB-TK classification for ${sample_base}"
echo "Filtered MAGs directory : $filtered_MAGs_dir"
echo "GTDB-TK output directory: $gtdbtk_output_dir"

# Set GTDB-TK reference database path
export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH}"

# Activate the GTDB-TK 2.7.1 conda environment
source activate /mnt/scratch2/igfs-anaconda/conda-envs/gtdbtk_2.7.1

mkdir -p "$gtdbtk_output_dir"

# ============================
# Run GTDB-TK classify_wf
# ============================
echo "Running GTDB-TK classify_wf..."

gtdbtk classify_wf \
  --genome_dir "$filtered_MAGs_dir" \
  --out_dir "$gtdbtk_output_dir" \
  --extension fa \
  --cpus 20 \
  --pplacer_cpus 20 \
  -x .fa


if [[ \$? -ne 0 ]]; then
  echo "ERROR: GTDB-TK classify_wf failed for ${sample_base}."
  source deactivate
  exit 1
fi

echo "GTDB-TK classify_wf completed for ${sample_base}."

# ============================
# Summarise results
# ============================
bac_summary="$gtdbtk_output_dir/classify/gtdbtk.bac120.summary.tsv"
ar_summary="$gtdbtk_output_dir/classify/gtdbtk.ar53.summary.tsv"

echo ""
echo "=========================================="
echo "GTDB-TK classification summary for ${sample_base}"
echo "  Input MAGs : \$(ls "$filtered_MAGs_dir"/*.fa | wc -l)"

if [[ -f "\$bac_summary" ]]; then
  bac_count=\$(tail -n +2 "\$bac_summary" | wc -l)
  echo "  Bacterial MAGs classified : \$bac_count"
fi

if [[ -f "\$ar_summary" ]]; then
  ar_count=\$(tail -n +2 "\$ar_summary" | wc -l)
  echo "  Archaeal MAGs classified  : \$ar_count"
fi

echo "  Results written to: $gtdbtk_output_dir"
echo "=========================================="

source deactivate

echo "Full GTDB-TK pipeline completed for ${sample_base}"

EOL

          echo "Submitting GTDB-TK job for $sample_base..."
          sbatch "$job_script"

        fi
      done
    fi
  done
done

echo "All GTDB-TK classification jobs submitted!"

