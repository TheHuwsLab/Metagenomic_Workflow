#!/bin/bash

#SBATCH --job-name=fastp_Ran
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --error=/users/3057556/jobs/fastp_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/fastp_Ran_%j.log
#SBATCH --partition=k2-hipri,k2-medpri,k2-lowpri,k2-bioinf
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Exit on errors, unset variables, and failed commands in pipelines
set -euo pipefail

# Load Anaconda
module load apps/anaconda3/2024.06/bin

# Activate the conda environment containing fastp
source activate /mnt/scratch2/igfs-anaconda/conda-envs/fastp_1.0.1

# -------------------------------------------------------------------------
# Working directory
# -------------------------------------------------------------------------

WORK_DIR="/mnt/scratch2/users/3057556/Ran/40-1310848928/00_fastq"

cd "$WORK_DIR" || {
    echo "ERROR: Could not change to working directory: $WORK_DIR"
    exit 1
}

echo "Working directory: $WORK_DIR"
echo "Starting fastp processing..."
echo

# -------------------------------------------------------------------------
# Sample directory prefixes to process
# -------------------------------------------------------------------------

DIR_TAGS=("3" "4" "P" "N")

# -------------------------------------------------------------------------
# Process each sample directory
# -------------------------------------------------------------------------

for TAG in "${DIR_TAGS[@]}"; do

    echo "============================================================"
    echo "Processing directories beginning with: $TAG"
    echo "============================================================"

    # Find directories beginning with the specified tag
    for dir in "$TAG"*; do

        # If no directory matches, skip
        [[ -d "$dir" ]] || continue

        echo
        echo "Processing directory: $dir"

        # -------------------------------------------------------------
        # Skip samples that have already completed successfully
        # -------------------------------------------------------------

        if [[ -f "$dir/fastp.succeed" ]]; then
            echo "Skipping $dir - fastp.succeed found"
            continue
        fi

        # -------------------------------------------------------------
        # Find R1 and R2 FASTQ files
        # -------------------------------------------------------------

        R1=$(find "$dir" -type f -name "*R1_001.fastq.gz" -print -quit)
        R2=$(find "$dir" -type f -name "*R2_001.fastq.gz" -print -quit)

        # Check that both files exist
        if [[ -z "$R1" || ! -f "$R1" ]]; then
            echo "WARNING: R1 file not found in $dir"
            continue
        fi

        if [[ -z "$R2" || ! -f "$R2" ]]; then
            echo "WARNING: R2 file not found in $dir"
            continue
        fi

        echo "R1: $R1"
        echo "R2: $R2"

        # -------------------------------------------------------------
        # Define output filenames
        # -------------------------------------------------------------

        out_R1="${R1/.fastq.gz/_trimmed_paired.fastq.gz}"
        out_R2="${R2/.fastq.gz/_trimmed_paired.fastq.gz}"

        out_unpaired_R1="${R1/.fastq.gz/_trimmed_unpaired.fastq.gz}"
        out_unpaired_R2="${R2/.fastq.gz/_trimmed_unpaired.fastq.gz}"

        report_html="${R1/.fastq.gz/_fastp_report.html}"

        report_json="${R1/.fastq.gz/_fastp_report.json}"

        echo "Output R1: $out_R1"
        echo "Output R2: $out_R2"
        echo

        # -------------------------------------------------------------
        # Run fastp
        # -------------------------------------------------------------

        if fastp \
            -i "$R1" \
            -I "$R2" \
            -o "$out_R1" \
            -O "$out_R2" \
            --unpaired1 "$out_unpaired_R1" \
            --unpaired2 "$out_unpaired_R2" \
            --html "$report_html" \
            --json "$report_json" \
            --detect_adapter_for_pe \
            --thread "$SLURM_CPUS_PER_TASK" \
            --qualified_quality_phred 20 \
            --length_required 50
        then

            # ---------------------------------------------------------
            # Mark sample as successfully completed
            # ---------------------------------------------------------

            touch "$dir/fastp.succeed"

            echo "SUCCESS: Finished processing $(basename "$dir")"

        else

            echo "ERROR: fastp failed for $dir"

        fi

    done

done

echo
echo "============================================================"
echo "All samples processed!"
echo "============================================================"