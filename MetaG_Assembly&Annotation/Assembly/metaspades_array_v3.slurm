#!/bin/bash

#SBATCH --job-name=Prepare_MetaSpades_Jobs_Ran
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --error=/users/3057556/jobs/Prepare_MetaSpades_Jobs_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_MetaSpades_Jobs_Ran_%j.log
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

JOB_DIR="/users/3057556/jobs/metaspades_jobs_Ran"

SPADES_ENV="/mnt/scratch2/igfs-anaconda/conda-envs/spades_4.2.0"

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
echo

# -------------------------------------------------------------------------
# Directory tags for Ran samples
# -------------------------------------------------------------------------

DIR_TAGS=("4" "N" "P")

# -------------------------------------------------------------------------
# Find samples and generate MetaSPAdes jobs
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
        # Find unmapped paired R1 FASTQ
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
        # Extract filenames
        # -------------------------------------------------------------

        filename=$(basename "$R1")

        # Example:
        #
        # SAMPLE_unmapped_trimmed_paired_R1.fastq.gz
        #
        # becomes:
        #
        # SAMPLE
        #
        sample_name="${filename%_unmapped_trimmed_paired_R1.fastq.gz}"

        # MetaSPAdes output directory
        outdir="${sample_name}_metaspades"

        # Full path to expected assembly
        contigs_file="$dir/$outdir/contigs.fasta"

        # Completion flag
        success_flag="$dir/${sample_name}_metaspades_assembly_done.flag"

        # Previous-step flag
        prev_flag="$dir/${sample_name}_decontamination_done.flag"

        # Job script
        job_script="$JOB_DIR/${sample_name}_metaspades_job.sh"

        # -------------------------------------------------------------
        # Print sample information
        # -------------------------------------------------------------

        echo "R1:            $R1"
        echo "R2:            $R2"
        echo "Sample:        $sample_name"
        echo "Output dir:    $outdir"
        echo "Contigs:       $contigs_file"
        echo "Previous flag: $prev_flag"
        echo "Success flag:  $success_flag"
        echo "Job script:    $job_script"

        # -------------------------------------------------------------
        # Check previous contamination step
        # -------------------------------------------------------------

        if [[ ! -f "$prev_flag" ]]; then
            echo "WARNING: Previous contamination step has not completed."
            echo "        Missing: $prev_flag"
            echo "        Skipping $sample_name."
            continue
        fi

        # -------------------------------------------------------------
        # Skip if MetaSPAdes has already completed
        # -------------------------------------------------------------

        if [[ -f "$success_flag" && -s "$contigs_file" ]]; then
            echo "Skipping $sample_name - MetaSPAdes already completed."
            continue
        fi

        # -------------------------------------------------------------
        # Also skip if contigs.fasta already exists and is non-empty
        # -------------------------------------------------------------

        if [[ -s "$contigs_file" ]]; then
            echo "Skipping $sample_name:"
            echo "  contigs.fasta already exists and is non-empty."
            continue
        fi

        # -------------------------------------------------------------
        # Generate SLURM MetaSPAdes job
        # -------------------------------------------------------------

        cat > "$job_script" <<EOL
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=300G
#SBATCH --job-name=MetaSpades_${sample_name}
#SBATCH --error=${JOB_DIR}/${sample_name}_%j.error
#SBATCH --output=${JOB_DIR}/${sample_name}_%j.log
#SBATCH --partition=k2-bioinf,k2-medpri,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# -------------------------------------------------------------------------
# Exit on errors, unset variables, and failed pipelines
# -------------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------------
# Load Anaconda and activate MetaSPAdes environment
# -------------------------------------------------------------------------

module load apps/anaconda3/2024.06/bin

source activate ${SPADES_ENV}

# -------------------------------------------------------------------------
# Define input/output files
# -------------------------------------------------------------------------

R1="${R1}"
R2="${R2}"

OUT_DIR="${dir}/${outdir}"

CONTIGS="\${OUT_DIR}/contigs.fasta"

PREV_FLAG="${prev_flag}"

SUCCESS_FLAG="${success_flag}"

# -------------------------------------------------------------------------
# Remove success flag if the job fails
# -------------------------------------------------------------------------

trap '
    rc=\$?
    if [[ \$rc -eq 0 ]]; then
        printf "%s\n" "\$(date --iso-8601=seconds)" > "\$SUCCESS_FLAG"
    else
        rm -f "\$SUCCESS_FLAG" 2>/dev/null || true
        echo "ERROR: MetaSPAdes failed for ${sample_name}" >&2
    fi
    exit \$rc
' EXIT

# -------------------------------------------------------------------------
# Verify previous contamination step completed
# -------------------------------------------------------------------------

if [[ ! -f "\$PREV_FLAG" ]]; then

    echo "ERROR: Previous-step flag not found:"
    echo "\$PREV_FLAG"
    echo "Cannot run MetaSPAdes."

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
# Run MetaSPAdes
# -------------------------------------------------------------------------

echo
echo "============================================================"
echo "Starting MetaSPAdes"
echo "============================================================"
echo "R1:       \$R1"
echo "R2:       \$R2"
echo "Output:   \$OUT_DIR"
echo "Threads:  \$SLURM_CPUS_PER_TASK"
echo "Memory:   300 GB"
echo

metaspades.py \
    -1 "\$R1" \
    -2 "\$R2" \
    -o "\$OUT_DIR" \
    -t "\$SLURM_CPUS_PER_TASK" \
    -m 300

# -------------------------------------------------------------------------
# Validate assembly
# -------------------------------------------------------------------------

if [[ ! -s "\$CONTIGS" ]]; then

    echo "ERROR: MetaSPAdes completed but contigs.fasta was not produced:"
    echo "\$CONTIGS"

    exit 1

fi

echo
echo "MetaSPAdes assembly completed successfully."
echo "Contigs:"
echo "\$CONTIGS"

echo
echo "Assembly size:"
du -h "\$CONTIGS"

echo
echo "Number of contigs:"
grep -c '^>' "\$CONTIGS" || true

echo
echo "MetaSPAdes job completed successfully."

EOL

        # -------------------------------------------------------------
        # Make generated script executable
        # -------------------------------------------------------------

        chmod +x "$job_script"

        # -------------------------------------------------------------
        # Submit job
        # -------------------------------------------------------------

        echo
        echo "Submitting MetaSPAdes job for $sample_name..."

        sbatch "$job_script"

    done

done

echo
echo "============================================================"
echo "All MetaSPAdes jobs submitted."
echo "============================================================"
