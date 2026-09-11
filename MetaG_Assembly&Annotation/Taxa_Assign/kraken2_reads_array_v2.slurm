#!/bin/bash

#SBATCH --job-name=Prepare_Kraken2_Reads_Jobs_Ran
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --error=/users/3057556/jobs/Prepare_Kraken2_Reads_Jobs_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_Kraken2_Reads_Jobs_Ran_%j.log
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# -------------------------------------------------------------------------
# Exit on errors, unset variables, and failed pipelines
# -------------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------------
# Load required modules
# -------------------------------------------------------------------------

module load apps/anaconda3/2024.06/bin

# -------------------------------------------------------------------------
# Define paths
# -------------------------------------------------------------------------

WORK_DIR="/mnt/scratch2/users/3057556/Ran/40-1310848928/00_fastq"

JOB_DIR="/users/3057556/jobs/Ran_kraken2_jobs"

KRAKEN_ENV="/mnt/scratch2/igfs-anaconda/conda-envs/kraken2_2.1.3"

DBNAME="/mnt/scratch2/igfs-databases/Kraken2_DBs/k2_pluspfp_20260626"

# -------------------------------------------------------------------------
# Create directory for generated SLURM jobs
# -------------------------------------------------------------------------

mkdir -p "$JOB_DIR"

# -------------------------------------------------------------------------
# Change to working directory
# -------------------------------------------------------------------------

cd "$WORK_DIR" || {
    echo "ERROR: Could not change to working directory:"
    echo "$WORK_DIR"
    exit 1
}

echo "Working directory: $WORK_DIR"
echo "Job directory:     $JOB_DIR"
echo "Kraken2 database:  $DBNAME"
echo

# -------------------------------------------------------------------------
# Check Kraken2 database
# -------------------------------------------------------------------------

if [[ ! -d "$DBNAME" ]]; then
    echo "ERROR: Kraken2 database does not exist:"
    echo "$DBNAME"
    exit 1
fi

# -------------------------------------------------------------------------
# Directory tags for Ran samples
# -------------------------------------------------------------------------

DIR_TAGS=("4" "N" "P")

# -------------------------------------------------------------------------
# Find samples and generate Kraken2 jobs
# -------------------------------------------------------------------------

for TAG in "${DIR_TAGS[@]}"; do

    echo
    echo "============================================================"
    echo "Searching directories beginning with: $TAG"
    echo "============================================================"

    for dir in "$TAG"*; do

        [[ -d "$dir" ]] || continue

        sample_dir=$(basename "$dir")

        echo
        echo "Processing directory: $dir"

        # -------------------------------------------------------------
        # Find decontaminated paired R1
        # -------------------------------------------------------------

        R1=$(find "$dir" \
            -maxdepth 1 \
            -type f \
            -name "*_unmapped_trimmed_paired_R1.fastq.gz" \
            -print -quit)

        if [[ -z "$R1" || ! -f "$R1" ]]; then
            echo "WARNING: No unmapped R1 FASTQ found in $dir"
            continue
        fi

        # -------------------------------------------------------------
        # Derive R2 from R1
        # -------------------------------------------------------------

        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

        if [[ ! -f "$R2" ]]; then
            echo "WARNING: R2 file not found:"
            echo "  Expected: $R2"
            continue
        fi

        # -------------------------------------------------------------
        # Extract sample name
        #
        # Example:
        #
        # 41644B_unmapped_trimmed_paired_R1.fastq.gz
        #
        # becomes:
        #
        # 41644B
        # -------------------------------------------------------------

        filename=$(basename "$R1")

        sample_name="${filename%_unmapped_trimmed_paired_R1.fastq.gz}"

        # -------------------------------------------------------------
        # Define Kraken2 output directory and files
        # -------------------------------------------------------------

        outdir="$dir/${sample_name}_kraken2"

        kraken_output="$outdir/${sample_name}_kraken2_reads.txt"

        report="$outdir/${sample_name}_kraken2_reads_report.txt"

        success_flag="$outdir/${sample_name}_kraken2_reads_done.flag"

        # -------------------------------------------------------------
        # Previous contamination-removal flag
        # -------------------------------------------------------------

        prev_flag="$dir/${sample_name}_decontamination_done.flag"

        # -------------------------------------------------------------
        # Generated job script
        # -------------------------------------------------------------

        job_script="$JOB_DIR/${sample_name}_kraken2_job.sh"

        # -------------------------------------------------------------
        # Print sample information
        # -------------------------------------------------------------

        echo "Sample:        $sample_name"
        echo "R1:            $R1"
        echo "R2:            $R2"
        echo "Output dir:    $outdir"
        echo "Kraken output: $kraken_output"
        echo "Report:        $report"
        echo "Previous flag: $prev_flag"
        echo "Success flag:  $success_flag"
        echo "Job script:    $job_script"

        # -------------------------------------------------------------
        # Check previous contamination step
        # -------------------------------------------------------------

        if [[ ! -f "$prev_flag" ]]; then
            echo
            echo "WARNING: Contamination removal has not completed."
            echo "Missing:"
            echo "  $prev_flag"
            echo "Skipping $sample_name."
            continue
        fi

        # -------------------------------------------------------------
        # Skip already completed Kraken2 samples
        # -------------------------------------------------------------

        if [[ -f "$success_flag" && -s "$report" ]]; then
            echo "Skipping $sample_name - Kraken2 already completed."
            continue
        fi

        # -------------------------------------------------------------
        # Generate the individual SLURM job
        # -------------------------------------------------------------

        cat > "$job_script" <<EOL
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=320G
#SBATCH --job-name=Kraken2_reads_${sample_name}
#SBATCH --error=${JOB_DIR}/${sample_name}_%j.error
#SBATCH --output=${JOB_DIR}/${sample_name}_%j.log
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# -------------------------------------------------------------------------
# Exit on errors, unset variables, and failed pipelines
# -------------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------------
# Load Anaconda and activate Kraken2 environment
# -------------------------------------------------------------------------

module load apps/anaconda3/2024.06/bin

source activate ${KRAKEN_ENV}

# -------------------------------------------------------------------------
# Define input/output paths
# -------------------------------------------------------------------------

R1="${R1}"
R2="${R2}"

OUT_DIR="${outdir}"

KRAKEN_OUTPUT="\${OUT_DIR}/${sample_name}_kraken2_reads.txt"

REPORT="\${OUT_DIR}/${sample_name}_kraken2_reads_report.txt"

PREV_FLAG="${prev_flag}"

SUCCESS_FLAG="${success_flag}"

DBNAME="${DBNAME}"

# -------------------------------------------------------------------------
# Remove success flag if the job fails
# -------------------------------------------------------------------------

trap '
    rc=\$?

    if [[ \$rc -eq 0 ]]; then

        printf "%s\n" "\$(date --iso-8601=seconds)" > "\$SUCCESS_FLAG"

    else

        rm -f "\$SUCCESS_FLAG" 2>/dev/null || true

        echo "ERROR: Kraken2 failed for ${sample_name}" >&2

    fi

    exit \$rc
' EXIT

# -------------------------------------------------------------------------
# Verify previous contamination-removal step
# -------------------------------------------------------------------------

if [[ ! -f "\$PREV_FLAG" ]]; then

    echo "ERROR: Previous-step flag not found:"
    echo "\$PREV_FLAG"

    exit 1

fi

echo "Previous-step flag found:"
echo "\$PREV_FLAG"

# -------------------------------------------------------------------------
# Verify input FASTQ files
# -------------------------------------------------------------------------

if [[ ! -s "\$R1" ]]; then

    echo "ERROR: R1 FASTQ not found or empty:"
    echo "\$R1"

    exit 1

fi

if [[ ! -s "\$R2" ]]; then

    echo "ERROR: R2 FASTQ not found or empty:"
    echo "\$R2"

    exit 1

fi

# -------------------------------------------------------------------------
# Create output directory
# -------------------------------------------------------------------------

mkdir -p "\$OUT_DIR"

# -------------------------------------------------------------------------
# Run Kraken2
# -------------------------------------------------------------------------

echo
echo "============================================================"
echo "Starting Kraken2"
echo "============================================================"

echo "Sample:   ${sample_name}"
echo "R1:       \$R1"
echo "R2:       \$R2"
echo "Database: \$DBNAME"
echo "Output:   \$OUT_DIR"
echo "Threads:  \$SLURM_CPUS_PER_TASK"
echo

kraken2 \
    --threads "\$SLURM_CPUS_PER_TASK" \
    --output "\$KRAKEN_OUTPUT" \
    --report "\$REPORT" \
    --db "\$DBNAME" \
    --use-names \
    --paired \
    "\$R1" \
    "\$R2"

# -------------------------------------------------------------------------
# Validate Kraken2 output
# -------------------------------------------------------------------------

if [[ ! -s "\$KRAKEN_OUTPUT" ]]; then

    echo "ERROR: Kraken2 read-level output was not created or is empty:"
    echo "\$KRAKEN_OUTPUT"

    exit 1

fi

if [[ ! -s "\$REPORT" ]]; then

    echo "ERROR: Kraken2 report was not created or is empty:"
    echo "\$REPORT"

    exit 1

fi

# -------------------------------------------------------------------------
# Report results
# -------------------------------------------------------------------------

echo
echo "Kraken2 completed successfully."

echo "Read-level output:"
echo "  \$KRAKEN_OUTPUT"

echo "Kraken2 report:"
echo "  \$REPORT"

echo
echo "Report summary:"
cat "\$REPORT"

echo
echo "============================================================"
echo "Kraken2 job completed successfully."
echo "============================================================"

EOL

        # -------------------------------------------------------------
        # Make generated job executable
        # -------------------------------------------------------------

        chmod +x "$job_script"

        # -------------------------------------------------------------
        # Submit individual job
        # -------------------------------------------------------------

        echo
        echo "Submitting Kraken2 job for $sample_name..."

        sbatch "$job_script"

    done

done

echo
echo "============================================================"
echo "All Kraken2 jobs submitted."
echo "============================================================"
