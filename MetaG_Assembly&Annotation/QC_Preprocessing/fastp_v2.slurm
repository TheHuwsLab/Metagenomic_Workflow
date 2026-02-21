#!/bin/bash

#SBATCH --job-name=fastp_
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --error=/users/3057556/jobs/fastp_PN0441_%j.error
#SBATCH --output=/users/3057556/jobs/fastp_PN0441_%j.log
#SBATCH --partition=k2-hipri,k2-medpri,k2-lowpri,k2-bioinf
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL


module load apps/anaconda3/2024.06/bin

# Activate the conda environment with fastp installed
source activate /mnt/scratch2/igfs-anaconda/conda-envs/fastp_1.0.1


# Define the parent directory containing subdirectories
WORK_DIR="/mnt/scratch2/users/3057556/TODO/PN0441/"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAGS=("E" "L" "P")

for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    if [[ -d "$dir" ]]; then
      echo "Processing directory: $dir"

      # Skip if already succeeded
      if [[ -f "$dir/fastp.succeed" ]]; then
        echo "Skipping $dir - fastp.succeed found"
        continue
      fi

      # Detect R1 and R2 files - May need to modify the name catches here
      R1=$(find "$dir" -type f -name "*R1_001.fastq.gz" | head -n 1)
      R2=$(find "$dir" -type f -name "*R2_001.fastq.gz" | head -n 1)

      if [[ -f "$R1" && -f "$R2" ]]; then
          # Define output filenames in the same directory
          out_R1="${R1/.fastq.gz/_trimmed_paired.fastq.gz}"
          out_R2="${R2/.fastq.gz/_trimmed_paired.fastq.gz}"
          out_unpaired_R1="${R1/.fastq.gz/_trimmed_unpaired.fastq.gz}"
          out_unpaired_R2="${R2/.fastq.gz/_trimmed_unpaired.fastq.gz}"
          report_html="${R1/.fastq.gz/_fastp_report.html}"

          echo "Read files to be trimmed: \"$R1\" & \"$R2\""

          # Run fastp with unpaired output enabled
          if fastp -i "$R1" -I "$R2" -o "$out_R1" -O "$out_R2" \
                --unpaired1 "$out_unpaired_R1" --unpaired2 "$out_unpaired_R2" \
                --html "$report_html" --detect_adapter_for_pe \
                --thread 20 --qualified_quality_phred 20 --length_required 50; then
              # Mark success so future runs skip this sample
              touch "$dir/fastp.succeed"
              echo "Finished processing: $(basename "$dir")"
          else
              echo "ERROR: fastp failed for $dir"
          fi
      else
          echo "WARNING: Missing R1 or R2 files in $dir, skipping..."
      fi
    fi
  done
done

echo "All samples processed!"





