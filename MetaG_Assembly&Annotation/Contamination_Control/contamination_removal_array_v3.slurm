#!/bin/bash

#SBATCH --job-name=Prepare_Ran_Contam_Jobs
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --error=/users/3057556/jobs/Prepare_Contam_Jobs_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_Contam_Jobs_Ran_%j.log
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Exit on errors, unset variables, and failed pipelines
set -euo pipefail

# -------------------------------------------------------------------------
# Load required modules
# -------------------------------------------------------------------------

module load apps/anaconda3/2024.06/bin

# -------------------------------------------------------------------------
# Define paths
# -------------------------------------------------------------------------

WORK_DIR="/mnt/scratch2/users/3057556/Ran/40-1310848928/00_fastq"

JOB_DIR="/users/3057556/jobs/Ran_contamination_jobs"

CONTAM_DB="/mnt/scratch2/igfs-databases/reference_genomes/rumen_contamination_genomes/rumen_contamination_bt2_db"

BOWTIE_ENV="/mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1"

# -------------------------------------------------------------------------
# Create job directory
# -------------------------------------------------------------------------

mkdir -p "$JOB_DIR"

# -------------------------------------------------------------------------
# Check working directory
# -------------------------------------------------------------------------

cd "$WORK_DIR" || {
    echo "ERROR: Could not change to working directory: $WORK_DIR"
    exit 1
}

# -------------------------------------------------------------------------
# Check contamination database exists
# -------------------------------------------------------------------------

if [[ ! -f "${CONTAM_DB}.1.bt2" && ! -f "${CONTAM_DB}.1.bt2l" ]]; then
    echo "ERROR: Bowtie2 contamination database not found:"
    echo "$CONTAM_DB"
    exit 1
fi

# -------------------------------------------------------------------------
# Directory tags to process
# -------------------------------------------------------------------------

DIR_TAGS=("4" "N" "P")

# -------------------------------------------------------------------------
# Find samples and generate contamination-check jobs
# -------------------------------------------------------------------------

for TAG in "${DIR_TAGS[@]}"; do

    echo
    echo "============================================================"
    echo "Searching directories beginning with: $TAG"
    echo "============================================================"

    for dir in "$TAG"*; do

        [[ -d "$dir" ]] || continue

        sample_name=$(basename "$dir")

        echo
        echo "Processing sample directory: $dir"

        # -------------------------------------------------------------
        # Find trimmed paired R1
        # -------------------------------------------------------------

        R1=$(find "$dir" \
            -maxdepth 1 \
            -type f \
            -name "*R1_001_trimmed_paired.fastq.gz" \
            -print -quit)

        if [[ -z "$R1" || ! -f "$R1" ]]; then
            echo "WARNING: No trimmed R1 found in $dir"
            continue
        fi

        # -------------------------------------------------------------
        # Derive R2 from R1
        # -------------------------------------------------------------

        R2="${R1/_R1_/_R2_}"

        if [[ ! -f "$R2" ]]; then
            echo "WARNING: R2 file not found:"
            echo "  Expected: $R2"
            continue
        fi

        filename=$(basename "$R1")

        # -------------------------------------------------------------
        # Derive stable sample base name
        #
        # Example:
        # SAMPLE_R1_001_trimmed_paired.fastq.gz
        #
        # becomes:
        # SAMPLE
        # -------------------------------------------------------------

        base="${filename%%_R1_*}"

        sorted_bam="${base}_sorted.bam"

        unmapped_R1="${base}_unmapped_trimmed_paired_R1.fastq.gz"
        unmapped_R2="${base}_unmapped_trimmed_paired_R2.fastq.gz"

        # -------------------------------------------------------------
        # Success flag
        # -------------------------------------------------------------

        success_flag="$dir/${sample_name}_decontamination_done.flag"

        # -------------------------------------------------------------
        # Job script
        # -------------------------------------------------------------

        job_script="$JOB_DIR/${sample_name}_contamination_check_job.sh"

        echo "R1:           $R1"
        echo "R2:           $R2"
        echo "Sample base:  $base"
        echo "Sorted BAM:   $sorted_bam"
        echo "Unmapped R1:  $unmapped_R1"
        echo "Unmapped R2:  $unmapped_R2"
        echo "Job script:   $job_script"
        echo "------------------------------------------------------------"

        # -------------------------------------------------------------
        # Do not submit if already completed
        # -------------------------------------------------------------

        if [[ -f "$success_flag" ]]; then
            echo "Skipping $sample_name - contamination check already completed"
            continue
        fi

        # -------------------------------------------------------------
        # Generate SLURM job
        # -------------------------------------------------------------

        cat > "$job_script" <<EOL
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=ContCheck_${sample_name}
#SBATCH --error=${JOB_DIR}/${sample_name}_%j.error
#SBATCH --output=${JOB_DIR}/${sample_name}_%j.log
#SBATCH --partition=k2-hipri,k2-bioinf
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# -------------------------------------------------------------------------
# Load modules and activate Bowtie2 environment
# -------------------------------------------------------------------------

module load apps/anaconda3/2024.06/bin

source activate ${BOWTIE_ENV}

# -------------------------------------------------------------------------
# Input/output paths
# -------------------------------------------------------------------------

R1="${R1}"
R2="${R2}"

OUT_DIR="${dir}"

SORTED_BAM="\${OUT_DIR}/${sorted_bam}"

UNMAPPED_R1="\${OUT_DIR}/${unmapped_R1}"
UNMAPPED_R2="\${OUT_DIR}/${unmapped_R2}"

FASTP_FLAG="\${OUT_DIR}/fastp.succeed"

SUCCESS_FLAG="${success_flag}"

CONTAM_DB="${CONTAM_DB}"

# -------------------------------------------------------------------------
# Remove success flag if this job exits unsuccessfully
# -------------------------------------------------------------------------

trap '
    rc=\$?
    if [[ \$rc -eq 0 ]]; then
        printf "%s\n" "\$(date --iso-8601=seconds)" > "\$SUCCESS_FLAG"
    else
        rm -f "\$SUCCESS_FLAG"
        echo "ERROR: Contamination-check job failed for ${sample_name}" >&2
    fi
    exit \$rc
' EXIT

# -------------------------------------------------------------------------
# Verify fastp completed successfully
# -------------------------------------------------------------------------

if [[ ! -f "\$FASTP_FLAG" ]]; then

    echo "ERROR: fastp success flag not found:"
    echo "\$FASTP_FLAG"
    echo "Aborting contamination check."

    exit 1

fi

echo "fastp success flag found:"
echo "\$FASTP_FLAG"

# -------------------------------------------------------------------------
# Verify input FASTQ files
# -------------------------------------------------------------------------

if [[ ! -f "\$R1" ]]; then
    echo "ERROR: R1 file not found:"
    echo "\$R1"
    exit 1
fi

if [[ ! -f "\$R2" ]]; then
    echo "ERROR: R2 file not found:"
    echo "\$R2"
    exit 1
fi

# -------------------------------------------------------------------------
# Run Bowtie2 against contamination database
#
# Reads that map to the contamination database are retained in the BAM.
# Unmapped reads are subsequently extracted and used for downstream work.
# -------------------------------------------------------------------------

echo
echo "Starting Bowtie2 contamination screening..."
echo "R1: \$R1"
echo "R2: \$R2"
echo "Database: \$CONTAM_DB"
echo

bowtie2 \
    --threads "\$SLURM_CPUS_PER_TASK" \
    -x "\$CONTAM_DB" \
    -1 "\$R1" \
    -2 "\$R2" \
    | samtools view -@ "\$SLURM_CPUS_PER_TASK" -b - \
    | samtools sort -@ "\$SLURM_CPUS_PER_TASK" -o "\$SORTED_BAM" -

# -------------------------------------------------------------------------
# Index BAM
# -------------------------------------------------------------------------

samtools index -@ "\$SLURM_CPUS_PER_TASK" -c "\$SORTED_BAM"

# -------------------------------------------------------------------------
# Extract reads where BOTH mates are unmapped
#
# -f 12:
#   read unmapped (0x4)
#   mate unmapped (0x8)
#
# -F 256:
#   exclude secondary alignments
#
# -F 2048:
#   exclude supplementary alignments
# -------------------------------------------------------------------------

echo
echo "Extracting paired unmapped reads..."

UNMAPPED_R1_TMP="\${UNMAPPED_R1%.gz}"
UNMAPPED_R2_TMP="\${UNMAPPED_R2%.gz}"

rm -f "\$UNMAPPED_R1_TMP" "\$UNMAPPED_R2_TMP"

samtools collate \
    -@ "$SLURM_CPUS_PER_TASK" \
    -u \
    -O \
    "$SORTED_BAM" \
| samtools fastq \
    -@ "$SLURM_CPUS_PER_TASK" \
    -f 12 \
    -F 0x900 \
    -1 "$UNMAPPED_R1_TMP" \
    -2 "$UNMAPPED_R2_TMP" \
    -0 /dev/null \
    -s /dev/null \
    -n



# -------------------------------------------------------------------------
# Compress unmapped FASTQ files
# -------------------------------------------------------------------------

echo
echo "Compressing unmapped FASTQ files..."

pigz \
    -p "\$SLURM_CPUS_PER_TASK" \
    -c "\$UNMAPPED_R1_TMP" > "\$UNMAPPED_R1"

rm -f "\$UNMAPPED_R1_TMP"

pigz \
    -p "\$SLURM_CPUS_PER_TASK" \
    -c "\$UNMAPPED_R2_TMP" > "\$UNMAPPED_R2"

rm -f "\$UNMAPPED_R2_TMP"

# -------------------------------------------------------------------------
# Basic output validation
# -------------------------------------------------------------------------


if [[ ! -s "\$SORTED_BAM" ]]; then
    echo "ERROR: Sorted BAM was not created or is empty."
    exit 1
fi

if [[ ! -f "\$UNMAPPED_R1" ]]; then
    echo "ERROR: Unmapped R1 FASTQ was not created."
    exit 1
fi

if [[ ! -f "\$UNMAPPED_R2" ]]; then
    echo "ERROR: Unmapped R2 FASTQ was not created."
    exit 1
fi

R1_COUNT=$(gzip -cd "$UNMAPPED_R1" | awk 'END {print NR/4}')
R2_COUNT=$(gzip -cd "$UNMAPPED_R2" | awk 'END {print NR/4}')

echo "Unmapped paired R1 reads: $R1_COUNT"
echo "Unmapped paired R2 reads: $R2_COUNT"

if [[ "$R1_COUNT" -ne "$R2_COUNT" ]]; then
    echo "ERROR: Unmapped R1/R2 files contain different numbers of reads." >&2
    exit 1
fi

# -------------------------------------------------------------------------
# Report statistics
# -------------------------------------------------------------------------

echo
echo "Contamination check completed."
echo "Sorted BAM:"
echo "  \$SORTED_BAM"

echo "Unmapped R1:"
echo "  \$UNMAPPED_R1"

echo "Unmapped R2:"
echo "  \$UNMAPPED_R2"

echo
echo "BAM statistics:"
samtools flagstat "\$SORTED_BAM"

echo
echo "Job completed successfully."

EOL

        # -------------------------------------------------------------
        # Make generated script executable
        # -------------------------------------------------------------

        chmod +x "$job_script"

        # -------------------------------------------------------------
        # Submit job
        # -------------------------------------------------------------

        echo "Submitting contamination-check job for $sample_name..."

        sbatch "$job_script"

    done

done

echo
echo "============================================================"
echo "All contamination-check jobs submitted."
echo "============================================================"
