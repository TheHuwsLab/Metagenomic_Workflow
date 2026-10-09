#!/bin/bash
#SBATCH --job-name=Prepare_ReadMap_Jobs_Ran
#SBATCH --error=/users/3057556/jobs/Prepare_ReadMap_Jobs_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_ReadMap_Jobs_Ran_%j.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
shopt -s nullglob

# =============================================================================
# Contig -> CDS read-mapping pipeline v3
# =============================================================================
# v3 is isolated from v2:
#   * all analysis outputs are written below OUTPUT_ROOT
#   * all generated child jobs/logs are written below JOB_DIR
#   * each sample gets a success flag only after every required output exists
#   * the final readmapping_v3.success flag is written only after every
#     submitted/previously-running sample reports success
#   * the child job is generated with a quoted heredoc so runtime variables
#     cannot be expanded by this preparation script
#   * rerunning the preparation script does not duplicate active jobs
#   * existing complete results self-heal their missing success flags
# =============================================================================

module load apps/anaconda3/2024.06/bin

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
WORK_DIR="/mnt/scratch2/users/3057556/Ran/40-1310848928/00_fastq"

JOB_DIR="/users/3057556/jobs/Ran_readmap_jobs"

OUTPUT_ROOT="${WORK_DIR}/readmapping"

PREP_SUCCESS_FLAG="${JOB_DIR}/Prepare_ReadMap.success"

GLOBAL_SUCCESS_FLAG="${OUTPUT_ROOT}/readmapping.success"

COMPLETION_JOB="${JOB_DIR}/ReadMap_completion_check.sh"

COMPLETION_LOG="${JOB_DIR}/ReadMap_completion_check.log"

COMPLETION_ERR="${JOB_DIR}/ReadMap_completion_check.error"


mkdir -p "$JOB_DIR" "$OUTPUT_ROOT"
cd "$WORK_DIR"

DIR_TAGS=("NC")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
log() {
    printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$*"
}

job_is_active() {
    local job_name="$1"
    command -v squeue >/dev/null 2>&1 || return 1
    squeue -h -u "$USER" -o '%j' 2>/dev/null | grep -Fxq "$job_name"
}

active_job_id() {
    local job_name="$1"
    command -v squeue >/dev/null 2>&1 || return 1
    squeue -h -u "$USER" -o '%i|%j' 2>/dev/null \
        | awk -F'|' -v name="$job_name" '$2 == name {print $1; exit}'
}

# -----------------------------------------------------------------------------
# Counters and state
# -----------------------------------------------------------------------------
submitted=0
skipped_complete=0
skipped_missing=0
skipped_active=0
failed_submission=0

# Arrays used to build the final completion-check dependency.
declared_samples=()
success_flags=()
dependency_job_ids=()

rm -f "$PREP_SUCCESS_FLAG" "$GLOBAL_SUCCESS_FLAG"

log "============================================================"
log "Contig -> CDS read-mapping pipeline v3"
log "============================================================"
log "Working directory: $WORK_DIR"
log "Output root:       $OUTPUT_ROOT"
log "Job directory:     $JOB_DIR"
log "Global flag:       $GLOBAL_SUCCESS_FLAG"
log "============================================================"

# -----------------------------------------------------------------------------
# Discover samples and prepare child jobs
# -----------------------------------------------------------------------------
for TAG in "${DIR_TAGS[@]}"; do
    for dir in "./${TAG}"*; do
        [[ -d "$dir" ]] || continue

        sample_name="$(basename "$dir")"

        read_files=("$dir"/*_unmapped_trimmed_R1.fastq*)
        if (( ${#read_files[@]} == 0 )); then
            log "Skipping $sample_name - no *_unmapped_trimmed_R1.fastq* file found."
            ((skipped_missing+=1))
            continue
        fi

        reads_file_1="${read_files[0]}"
        reads_file_2="${reads_file_1/_R1/_R2}"

        contigs="${dir}/${sample_name}_unmapped_trimmed_metaspades/contigs.fasta"
        gff="${dir}/${sample_name}_pyrodigal/${sample_name}_pyrodigal.gff"

        # Every v3 output for this sample lives here. v2 output is untouched.
        outdir="${OUTPUT_ROOT}/${sample_name}"
        bowtie_db="${outdir}/${sample_name}_bowtie2_db"
        bam="${outdir}/${sample_name}_readmapped_sorted.bam"
        bai="${bam}.bai"
        contig_summary="${outdir}/${sample_name}_readmapped_contig_summary.txt"
        cds_summary="${outdir}/${sample_name}_readmapped_cds_summary.txt"
        success_flag="${outdir}/${sample_name}_readmapping.success"

        job_name="ReadMap_v3_${sample_name}"
        job_script="${JOB_DIR}/${sample_name}_readmap_v3_job.sh"
        job_stdout="${JOB_DIR}/${sample_name}_output.log"
        job_stderr="${JOB_DIR}/${sample_name}_error.log"

        log ""
        log "############################################################"
        log "Sample: $sample_name"
        log "############################################################"
        log "R1:      $reads_file_1"
        log "R2:      $reads_file_2"
        log "Contigs: $contigs"
        log "GFF:     $gff"
        log "Output:  $outdir"

        # Record this sample only if all required inputs are available.
        if [[ ! -s "$reads_file_2" ]]; then
            log "Skipping $sample_name - paired R2 file not found: $reads_file_2"
            ((skipped_missing+=1))
            continue
        fi

        if [[ ! -s "$contigs" ]]; then
            log "Skipping $sample_name - contigs not found: $contigs"
            ((skipped_missing+=1))
            continue
        fi

        if [[ ! -s "$gff" ]]; then
            log "Skipping $sample_name - Pyrodigal GFF not found: $gff"
            ((skipped_missing+=1))
            continue
        fi

        declared_samples+=("$sample_name")
        success_flags+=("$success_flag")

        # Fully complete sample: do not touch it.
        if [[ -s "$bam" && -s "$bai" && -s "$contig_summary" && -s "$cds_summary" && -f "$success_flag" ]]; then
            log "$sample_name already completed successfully."
            log "  Success flag: $success_flag"
            ((skipped_complete+=1))
            continue
        fi

        # Self-heal a missing flag when every result file is already present.
        if [[ -s "$bam" && -s "$bai" && -s "$contig_summary" && -s "$cds_summary" && ! -f "$success_flag" ]]; then
            log "All v3 outputs already exist for $sample_name; creating missing success flag."
            mkdir -p "$outdir"
            printf '%s\n' "$(date --iso-8601=seconds)" > "$success_flag"
            ((skipped_complete+=1))
            continue
        fi

        # Do not duplicate an active child job if preparation is rerun.
        if job_is_active "$job_name"; then
            existing_job_id="$(active_job_id "$job_name" || true)"
            log "Skipping submission for $sample_name - active job found: ${job_name} (ID ${existing_job_id:-unknown})"
            if [[ -n "${existing_job_id:-}" ]]; then
                dependency_job_ids+=("$existing_job_id")
            fi
            ((skipped_active+=1))
            continue
        fi

        rm -f "$success_flag"
        mkdir -p "$outdir"

        # ---------------------------------------------------------------------
        # Generate the child job.
        # Static values are emitted as shell assignments. The runtime body is
        # deliberately protected by a quoted heredoc.
        # ---------------------------------------------------------------------
        {
            cat <<HEADER
#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --job-name=${job_name}
#SBATCH --error=${job_stderr}
#SBATCH --output=${job_stdout}
#SBATCH --partition=k2-medpri,k2-bioinf,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

WORK_DIR=$(printf '%q' "$WORK_DIR")
SAMPLE_NAME=$(printf '%q' "$sample_name")
READS_1=$(printf '%q' "$reads_file_1")
READS_2=$(printf '%q' "$reads_file_2")
CONTIGS=$(printf '%q' "$contigs")
GFF=$(printf '%q' "$gff")
OUTDIR=$(printf '%q' "$outdir")
BOWTIE_DB=$(printf '%q' "$bowtie_db")
BAM=$(printf '%q' "$bam")
BAI=$(printf '%q' "$bai")
CONTIG_SUMMARY=$(printf '%q' "$contig_summary")
CDS_SUMMARY=$(printf '%q' "$cds_summary")
SUCCESS_FLAG=$(printf '%q' "$success_flag")

HEADER
            cat <<'JOB_BODY'
cleanup_on_failure() {
    local exit_code="${1:-1}"

    # Do not allow the EXIT trap to fire again while exiting cleanup.
    trap - EXIT

    echo >&2
    echo "============================================================" >&2
    echo "ERROR: Read-mapping pipeline failed for ${SAMPLE_NAME}" >&2
    echo "Exit code: ${exit_code}" >&2
    echo "Output directory: ${OUTDIR}" >&2
    echo "============================================================" >&2

    rm -f "$SUCCESS_FLAG"
    exit "$exit_code"
}

on_exit() {
    local exit_code="${1:-0}"

    if [[ "$exit_code" -ne 0 ]]; then
        cleanup_on_failure "$exit_code"
    fi
}

# Capture the status before entering the function. This is safe with set -u.
trap 'rc=$?; on_exit "$rc"' EXIT

source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

mkdir -p "$OUTDIR"

# Build Bowtie2 database if it does not already exist.
if [[ ! -f "${BOWTIE_DB}.1.bt2" && ! -f "${BOWTIE_DB}.1.bt2l" ]]; then
    echo "Building Bowtie2 database..."
    bowtie2-build "$CONTIGS" "$BOWTIE_DB"
fi

if [[ ! -f "${BOWTIE_DB}.1.bt2" && ! -f "${BOWTIE_DB}.1.bt2l" ]]; then
    echo "ERROR: Bowtie2 index was not created." >&2
    exit 1
fi

echo "Mapping reads..."
bowtie2 \
    -x "$BOWTIE_DB" \
    -1 "$READS_1" \
    -2 "$READS_2" \
    --very-sensitive-local \
    -p "${SLURM_CPUS_PER_TASK}" \
    | samtools view -bS -@ "${SLURM_CPUS_PER_TASK}" - \
    | samtools sort -@ "${SLURM_CPUS_PER_TASK}" \
        -T "${OUTDIR}/${SAMPLE_NAME}_sorttmp" \
        -o "$BAM"

[[ -s "$BAM" ]] || { echo "ERROR: BAM was not created." >&2; exit 1; }

echo "Indexing BAM..."
samtools index -@ "${SLURM_CPUS_PER_TASK}" "$BAM"
[[ -s "$BAI" ]] || { echo "ERROR: BAM index was not created." >&2; exit 1; }

echo "Generating contig summary..."
total_reads="$(samtools view -c "$BAM")"
samtools idxstats "$BAM" \
    | awk -v total="$total_reads" '
BEGIN {
    print "Contig_Name\tContig_Length\tMapped_Reads\tTotal_Reads\tProportion"
}
{
    mapped = $3
    len = $2
    if (total > 0) prop = mapped / total
    else prop = 0
    print $1 "\t" len "\t" mapped "\t" total "\t" prop
}' > "$CONTIG_SUMMARY"

[[ -s "$CONTIG_SUMMARY" ]] || { echo "ERROR: Contig summary was not created." >&2; exit 1; }

echo "Generating CDS summary..."
bedtools intersect \
    -abam "$BAM" \
    -b "$GFF" \
    -wa \
    -wb \
    -bed \
    | awk '
{
    attr = $NF
    if (match(attr, /ID=([^;]+)/, m)) print m[1]
    else if (match(attr, /gene=([^;]+)/, m)) print m[1]
    else if (match(attr, /Parent=([^;]+)/, m)) print m[1]
    else if (match(attr, /Name=([^;]+)/, m)) print m[1]
    else {
        sub(/;.*/, "", attr)
        print attr
    }
}' \
    | sort \
    | uniq -c > "$CDS_SUMMARY"

[[ -s "$CDS_SUMMARY" ]] || { echo "ERROR: CDS summary was not created." >&2; exit 1; }

rm -f "${OUTDIR}/${SAMPLE_NAME}_sorttmp"* 2>/dev/null || true

# Write the sample success flag only after every required result exists.
printf '%s\n' "$(date --iso-8601=seconds)" > "$SUCCESS_FLAG"

echo
echo "============================================================"
echo "Read-mapping completed successfully for ${SAMPLE_NAME}"
echo "Success flag: ${SUCCESS_FLAG}"
echo "============================================================"
JOB_BODY
        } > "$job_script"

        chmod 750 "$job_script"

        if ! bash -n "$job_script"; then
            log "ERROR: Generated child job failed bash -n: $job_script"
            rm -f "$job_script"
            rm -f "$success_flag"
            ((failed_submission+=1))
            continue
        fi

        log "Submitting job for $sample_name..."
        if sbatch_output=$(sbatch "$job_script" 2>&1); then
            log "$sbatch_output"
            child_job_id="$(sed -n 's/^Submitted batch job \([0-9][0-9]*\).*$/\1/p' <<< "$sbatch_output" | head -n1)"
            if [[ -n "$child_job_id" ]]; then
                dependency_job_ids+=("$child_job_id")
            fi
            ((submitted+=1))
        else
            rc=$?
            log "ERROR: sbatch submission failed for $sample_name (exit $rc)"
            log "$sbatch_output"
            rm -f "$success_flag"
            ((failed_submission+=1))
        fi
    done
done

# -----------------------------------------------------------------------------
# Preparation summary
# -----------------------------------------------------------------------------
log ""
log "============================================================"
log "Preparation summary"
log "============================================================"
log "Candidate samples:   ${#declared_samples[@]}"
log "Submitted:            $submitted"
log "Already complete:     $skipped_complete"
log "Already active:       $skipped_active"
log "Missing inputs:       $skipped_missing"
log "Submission failures:  $failed_submission"
log "============================================================"

if (( failed_submission > 0 )); then
    log "Preparation failed; no preparation success flag written."
    exit 1
fi

# Preparation success means all requested child jobs were either successfully
# submitted, already complete, or already active. It does not mean the analysis
# has finished.
printf '%s\n' "$(date --iso-8601=seconds)" > "$PREP_SUCCESS_FLAG"
log "Preparation success flag: $PREP_SUCCESS_FLAG"

# -----------------------------------------------------------------------------
# Build and submit a lightweight final completion checker.
# It has no heavy resources and writes GLOBAL_SUCCESS_FLAG only after all
# candidate samples have their sample-level .success flags.
# -----------------------------------------------------------------------------
{
    cat <<CHECK_HEADER
#!/bin/bash
#SBATCH --job-name=ReadMap_v3_FinalCheck
#SBATCH --error=${COMPLETION_ERR}
#SBATCH --output=${COMPLETION_LOG}
#SBATCH --partition=k2-hipri
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail

OUTPUT_ROOT=$(printf '%q' "$OUTPUT_ROOT")
GLOBAL_SUCCESS_FLAG=$(printf '%q' "$GLOBAL_SUCCESS_FLAG")
CHECK_HEADER
    printf 'EXPECTED_FLAGS=('
    printf '%q ' "${success_flags[@]}"
    printf ')\n'
    cat <<'CHECK_BODY'

rm -f "$GLOBAL_SUCCESS_FLAG"

missing=0
for flag in "${EXPECTED_FLAGS[@]}"; do
    if [[ ! -f "$flag" ]]; then
        echo "Missing sample success flag: $flag" >&2
        missing=1
    fi
done

if [[ "$missing" -ne 0 ]]; then
    echo "Read-mapping v3 did not complete successfully for all samples." >&2
    exit 1
fi

printf '%s\n' "$(date --iso-8601=seconds)" > "$GLOBAL_SUCCESS_FLAG"
echo "============================================================"
echo "READ-MAPPING V3 COMPLETE"
echo "Global success flag: $GLOBAL_SUCCESS_FLAG"
echo "============================================================"
CHECK_BODY
} > "$COMPLETION_JOB"

chmod 750 "$COMPLETION_JOB"

if ! bash -n "$COMPLETION_JOB"; then
    log "ERROR: Completion checker failed bash -n: $COMPLETION_JOB"
    exit 1
fi

# If there is nothing left to wait for, the checker can run immediately.
# Otherwise make it depend on every submitted/known-active child job.
if ((${#dependency_job_ids[@]} > 0)); then
    dependency_list="$(IFS=:; echo "${dependency_job_ids[*]}")"
    log "Submitting final completion checker with dependency: $dependency_list"
    sbatch --dependency="afterok:${dependency_list}" "$COMPLETION_JOB"
else
    log "No active/submitted child jobs remain; running completion checker immediately."
    if ! sbatch "$COMPLETION_JOB"; then
        log "ERROR: Failed to submit final completion checker."
        exit 1
    fi
fi