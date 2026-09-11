#!/bin/bash

#SBATCH --job-name=Prepare_Pyrodigal_Jobs_Ran
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --error=/users/3057556/jobs/Prepare_Pyrodigal_Jobs_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_Pyrodigal_Jobs_Ran_%j.log
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

JOB_DIR="/users/3057556/jobs/Ran_pyrodigal_jobs"

PYRODIGAL_ENV="/mnt/scratch2/igfs-anaconda/conda-envs/pyrodigal_3.6.3"

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
# Find MetaSPAdes assemblies and generate Pyrodigal jobs
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
        # Find MetaSPAdes assembly
        # -------------------------------------------------------------

        contigs_file=$(find "$dir" \
            -maxdepth 2 \
            -type f \
            -path "*/_metaspades/contigs.fasta" \
            -print -quit)

        # The above pattern may not match sample-specific directory names,
        # so use a broader search if necessary.
        if [[ -z "$contigs_file" ]]; then
            contigs_file=$(find "$dir" \
                -maxdepth 2 \
                -type f \
                -name "contigs.fasta" \
                -path "*_metaspades/*" \
                -print -quit)
        fi

        if [[ -z "$contigs_file" || ! -s "$contigs_file" ]]; then
            echo "WARNING: No non-empty MetaSPAdes contigs.fasta found in $dir"
            continue
        fi

        # -------------------------------------------------------------
        # Determine sample name
        #
        # Example:
        #
        # 41644B/41644B_metaspades/contigs.fasta
        #
        # becomes:
        #
        # 41644B
        # -------------------------------------------------------------

        sample_name=$(basename "$(dirname "$contigs_file")")
        sample_name="${sample_name%_metaspades}"

        # -------------------------------------------------------------
        # Define Pyrodigal output directory/files
        # -------------------------------------------------------------

        outdir="$dir/${sample_name}_pyrodigal"

        pyrodigal_out_gff="$outdir/${sample_name}_pyrodigal.gff"
        pyrodigal_out_nt="$outdir/${sample_name}_pyrodigal_nt.fa"
        pyrodigal_out_aa="$outdir/${sample_name}_pyrodigal_aa.fa"

        # -------------------------------------------------------------
        # Define success flag
        # -------------------------------------------------------------

        success_flag="$outdir/${sample_name}_pyrodigal_done.flag"

        # -------------------------------------------------------------
        # Define MetaSPAdes completion flag
        # -------------------------------------------------------------

        metaspades_flag="$dir/${sample_name}_metaspades_assembly_done.flag"

        # -------------------------------------------------------------
        # Define generated job script
        # -------------------------------------------------------------

        job_script="$JOB_DIR/${sample_name}_pyrodigal_job.sh"

        # -------------------------------------------------------------
        # Print sample information
        # -------------------------------------------------------------

        echo "Sample:          $sample_name"
        echo "Assembly:        $contigs_file"
        echo "Output directory: $outdir"
        echo "GFF:             $pyrodigal_out_gff"
        echo "Nucleotide:      $pyrodigal_out_nt"
        echo "Protein:         $pyrodigal_out_aa"
        echo "MetaSPAdes flag: $metaspades_flag"
        echo "Success flag:    $success_flag"
        echo "Job script:      $job_script"

        # -------------------------------------------------------------
        # Require successful MetaSPAdes completion
        # -------------------------------------------------------------

        if [[ ! -f "$metaspades_flag" ]]; then
            echo
            echo "WARNING: MetaSPAdes completion flag not found:"
            echo "  $metaspades_flag"
            echo "Skipping $sample_name."
            continue
        fi

        # -------------------------------------------------------------
        # Skip if Pyrodigal already completed
        # -------------------------------------------------------------

        if [[ -f "$success_flag" \
              && -s "$pyrodigal_out_gff" \
              && -s "$pyrodigal_out_nt" \
              && -s "$pyrodigal_out_aa" ]]; then

            echo "Skipping $sample_name - Pyrodigal already completed."
            continue

        fi

        # -------------------------------------------------------------
        # Skip if all outputs exist even if the flag is missing
        # -------------------------------------------------------------

        if [[ -s "$pyrodigal_out_gff" \
              && -s "$pyrodigal_out_nt" \
              && -s "$pyrodigal_out_aa" ]]; then

            echo "All Pyrodigal outputs already exist for $sample_name."
            echo "Skipping job generation."

            continue

        fi

        # -------------------------------------------------------------
        # Generate individual SLURM job
        # -------------------------------------------------------------

        cat > "$job_script" <<EOL
#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --job-name=Pyrodigal_${sample_name}
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
# Load Anaconda and activate Pyrodigal environment
# -------------------------------------------------------------------------

module load apps/anaconda3/2024.06/bin

source activate ${PYRODIGAL_ENV}

# -------------------------------------------------------------------------
# Define paths
# -------------------------------------------------------------------------

CONTIGS="${contigs_file}"

OUT_DIR="${outdir}"

GFF="\${OUT_DIR}/${sample_name}_pyrodigal.gff"

NT="\${OUT_DIR}/${sample_name}_pyrodigal_nt.fa"

AA="\${OUT_DIR}/${sample_name}_pyrodigal_aa.fa"

PREV_FLAG="${metaspades_flag}"

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

        echo "ERROR: Pyrodigal failed for ${sample_name}" >&2

    fi

    exit \$rc
' EXIT

# -------------------------------------------------------------------------
# Verify MetaSPAdes completed successfully
# -------------------------------------------------------------------------

if [[ ! -f "\$PREV_FLAG" ]]; then

    echo "ERROR: MetaSPAdes completion flag not found:"
    echo "\$PREV_FLAG"

    exit 1

fi

echo "MetaSPAdes completion flag found:"
echo "\$PREV_FLAG"

# -------------------------------------------------------------------------
# Verify input assembly
# -------------------------------------------------------------------------

if [[ ! -s "\$CONTIGS" ]]; then

    echo "ERROR: contigs.fasta not found or empty:"
    echo "\$CONTIGS"

    exit 1

fi

# -------------------------------------------------------------------------
# Create output directory
# -------------------------------------------------------------------------

mkdir -p "\$OUT_DIR"

# -------------------------------------------------------------------------
# Run Pyrodigal
# -------------------------------------------------------------------------

echo
echo "============================================================"
echo "Starting Pyrodigal"
echo "============================================================"

echo "Sample:      ${sample_name}"
echo "Input:       \$CONTIGS"
echo "Output dir:  \$OUT_DIR"
echo "Threads:     \$SLURM_CPUS_PER_TASK"
echo "Mode:        meta"

echo

pyrodigal \
    -i "\$CONTIGS" \
    -o "\$GFF" \
    -a "\$AA" \
    -d "\$NT" \
    -p meta \
    -j "\$SLURM_CPUS_PER_TASK"

# -------------------------------------------------------------------------
# Validate outputs
# -------------------------------------------------------------------------

if [[ ! -s "\$GFF" ]]; then

    echo "ERROR: Pyrodigal GFF output was not created or is empty:"
    echo "\$GFF"

    exit 1

fi

if [[ ! -s "\$NT" ]]; then

    echo "ERROR: Pyrodigal nucleotide output was not created or is empty:"
    echo "\$NT"

    exit 1

fi

if [[ ! -s "\$AA" ]]; then

    echo "ERROR: Pyrodigal protein output was not created or is empty:"
    echo "\$AA"

    exit 1

fi

# -------------------------------------------------------------------------
# Basic output statistics
# -------------------------------------------------------------------------

echo
echo "Pyrodigal completed successfully."

echo
echo "Output files:"
echo "  GFF: \$GFF"
echo "  NT:  \$NT"
echo "  AA:  \$AA"

echo
echo "Protein sequences:"
grep -c '^>' "\$AA" || true

echo
echo "Nucleotide sequences:"
grep -c '^>' "\$NT" || true

echo
echo "============================================================"
echo "Pyrodigal job completed successfully."
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
        echo "Submitting Pyrodigal job for $sample_name..."

        sbatch "$job_script"

    done

done

echo
echo "============================================================"
echo "All Pyrodigal jobs submitted."
echo "============================================================"
