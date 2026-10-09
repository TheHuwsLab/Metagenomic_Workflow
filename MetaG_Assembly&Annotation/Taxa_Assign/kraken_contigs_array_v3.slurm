#!/bin/bash
#SBATCH --job-name=Prepare_Kraken2_Contigs_Ran
#SBATCH --error=/users/3057556/jobs/Prepare_Kraken2_Contigs_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_Kraken2_Contigs_Ran_%j.log
#SBATCH --partition=k2-hipri
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail
shopt -s nullglob

# =============================================================================
# Kraken2 contig classification + aggregation workflow for Ran.
#
# Individual results stay in each sample's <sample>_kraken2 directory.
# Combined results are written to:
#   $WORK_DIR/Kraken2_contigs_combined/
#
# The script is restartable:
#   - Current, complete sample outputs and flags are reused.
#   - Missing sample flags are healed from current, complete outputs.
#   - Existing completed classifications are not rerun just to combine them.
#   - Incomplete samples are submitted/retried, unless an equivalent job is active.
#   - The final checker waits for all submitted/active jobs (afterany), then
#     aggregates successful samples and records failed samples as warnings.
#   - A clean success flag is only written if all discovered samples succeeded.
#   - A completed-with-warnings flag is written if one or more samples failed.
#
# The aggregation is embedded in the generated checker; no separate Python file
# is needed. The only external per-sample commands are Kraken2 and seqtk.
# =============================================================================

module load apps/anaconda3/2024.06/bin

# Paths - Ran dataset
WORK_DIR="/mnt/scratch2/users/3057556/Ran/40-1310848928/00_fastq"
JOB_DIR="/users/3057556/jobs/Ran_kraken2_contig_jobs"
DBNAME="/mnt/scratch2/igfs-databases/Kraken2_DBs/k2_pluspfp_20260626"
KRAKEN_ENV="/mnt/scratch2/igfs-anaconda/conda-envs/kraken2_2.17.1"
DIR_TAGS=("4" "N" "P" "3")
SAMPLE_FLAG_SUFFIX="_kraken2.success"

# Combined output layout mirrors the read-level Kraken2 combiner.
OUTPUT_DIR="${WORK_DIR}/Kraken2_contigs_combined"
COMBINED_DIR="${OUTPUT_DIR}/combined"
SUMMARY_DIR="${OUTPUT_DIR}/summaries"
ABUNDANCE_DIR="${OUTPUT_DIR}/abundance"
VIS_DIR="${OUTPUT_DIR}/visualisations"
KRONA_DIR="${OUTPUT_DIR}/Krona"

PREP_SUCCESS_FLAG="${JOB_DIR}/Prepare_Kraken2_Contigs.success"
GLOBAL_SUCCESS_FLAG="${WORK_DIR}/Kraken2_contigs.success"
COMBINE_SUCCESS_FLAG="${OUTPUT_DIR}/Kraken2_contigs_combination.success"
COMBINE_WARNING_FLAG="${OUTPUT_DIR}/completed_with_warnings.flag"
SAMPLE_MANIFEST="${JOB_DIR}/Kraken2_contigs_samples.tsv"
COMPLETION_JOB="${JOB_DIR}/Kraken2_Contigs_completion_check.sh"
COMPLETION_LOG="${JOB_DIR}/Kraken2_Contigs_completion_check.log"
COMPLETION_ERR="${JOB_DIR}/Kraken2_Contigs_completion_check.error"

mkdir -p "$JOB_DIR" "$OUTPUT_DIR" "$COMBINED_DIR" "$SUMMARY_DIR" "$ABUNDANCE_DIR" "$VIS_DIR" "$KRONA_DIR"
cd "$WORK_DIR"

log() { printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$*"; }

job_is_active() {
    local sample_name="$1"
    command -v squeue >/dev/null 2>&1 || return 1
    squeue -h -u "${USER:-}" -o '%j' 2>/dev/null \
      | awk -v s="$sample_name" '$0 == "Kraken2_" s || $0 == "Kraken2_" s {found=1} END {exit !found}'
}

active_job_id() {
    local sample_name="$1"
    command -v squeue >/dev/null 2>&1 || return 1
    squeue -h -u "${USER:-}" -o '%i|%j' 2>/dev/null \
      | awk -F'|' -v s="$sample_name" '$2 == "Kraken2_" s || $2 == "Kraken2_" s {print $1; exit}'
}

jobid_file() { printf '%s/%s.kraken2.jobid\n' "$JOB_DIR" "$1"; }

outputs_complete() {
    local sample_name="$1" sample_dir="$2" outdir
    outdir="${sample_dir}/${sample_name}_kraken2"
    [[ -s "${outdir}/${sample_name}_kraken2.txt" && \
       -s "${outdir}/${sample_name}_kraken2_report.txt" && \
       -s "${outdir}/${sample_name}_kraken2_report_mpa.txt" ]]
}

outputs_current() {
    local sample_name="$1" sample_dir="$2" contigs="$3" outdir
    outdir="${sample_dir}/${sample_name}_kraken2"
    outputs_complete "$sample_name" "$sample_dir" || return 1
    [[ "$contigs" -ot "${outdir}/${sample_name}_kraken2.txt" ]] || return 1
    [[ "$contigs" -ot "${outdir}/${sample_name}_kraken2_report.txt" ]] || return 1
    [[ "$contigs" -ot "${outdir}/${sample_name}_kraken2_report_mpa.txt" ]] || return 1
    [[ "${outdir}/${sample_name}_kraken2_report.txt" -ot "${outdir}/${sample_name}_kraken2_report_mpa.txt" ]] || return 1
}

heal_sample_flag() {
    local sample_name="$1" sample_dir="$2" flag
    flag="${sample_dir}/${sample_name}${SAMPLE_FLAG_SUFFIX}"
    if [[ ! -f "$flag" ]]; then
        printf '%s\n' "$(date --iso-8601=seconds)" > "$flag"
        log "[HEAL] Recreated sample flag: $flag"
    else
        log "[OK] Sample flag already exists: $flag"
    fi
}

cleanup_on_failure() {
    local rc="${1:-1}"
    trap - EXIT
    if [[ "$rc" -ne 0 ]]; then
        rm -f "$PREP_SUCCESS_FLAG"
        log "ERROR: Kraken2 preparation failed (exit $rc)."
    fi
    exit "$rc"
}
trap 'cleanup_on_failure "$?"' EXIT

submitted=0
skipped_complete=0
skipped_active=0
healed_flags=0
failed_submission=0
missing_inputs=0
declared_samples=0
dependency_job_ids=()

log "============================================================"
log "Preparing Kraken2 contig jobs and combined summaries"
log "Working directory: $WORK_DIR"
log "Job directory: $JOB_DIR"
log "Combined output directory: $OUTPUT_DIR"
log "Database: $DBNAME"
log "============================================================"

# Build the complete manifest first. Samples without contigs remain in the
# manifest so the final combined run can report them as warnings.
printf 'Sample\tSample_directory\tContig_file\tSuccess_flag\n' > "${SAMPLE_MANIFEST}.tmp"
for TAG in "${DIR_TAGS[@]}"; do
    for dir in "./${TAG}"*; do
        [[ -d "$dir" ]] || continue
        sample_name="$(basename "$dir")"
        sample_dir="${WORK_DIR}/${dir#./}"
        contigs="${sample_dir}/${sample_name}_metaspades/contigs.fasta"
        sample_flag="${sample_dir}/${sample_name}${SAMPLE_FLAG_SUFFIX}"
        printf '%s\t%s\t%s\t%s\n' "$sample_name" "$sample_dir" "$contigs" "$sample_flag" >> "${SAMPLE_MANIFEST}.tmp"
        ((declared_samples+=1))
    done
done
mv "${SAMPLE_MANIFEST}.tmp" "$SAMPLE_MANIFEST"

if (( declared_samples == 0 )); then
    log "ERROR: No sample directories found for tags: ${DIR_TAGS[*]}"
    exit 1
fi

# Prepare/submit only samples that are not already current and complete.
while IFS=$'\t' read -r sample_name sample_dir contigs success_flag; do
    [[ "$sample_name" == "Sample" || -z "$sample_name" ]] && continue
    outdir="${sample_dir}/${sample_name}_kraken2"
    jid_file="$(jobid_file "$sample_name")"
    job_name="Kraken2_${sample_name}"
    job_script="${JOB_DIR}/${sample_name}_kraken2_contig_job.sh"
    job_stdout="${JOB_DIR}/${sample_name}_output.log"
    job_stderr="${JOB_DIR}/${sample_name}_error.log"

    if [[ ! -s "$contigs" ]]; then
        log "[WARNING] $sample_name has no contig FASTA: $contigs"
        ((missing_inputs+=1))
        continue
    fi

    mkdir -p "$outdir"

    if outputs_current "$sample_name" "$sample_dir" "$contigs" && [[ -f "$success_flag" ]]; then
        log "[SKIP] $sample_name already has current outputs and its success flag."
        ((skipped_complete+=1))
        continue
    fi

    if outputs_current "$sample_name" "$sample_dir" "$contigs" && [[ ! -f "$success_flag" ]]; then
        heal_sample_flag "$sample_name" "$sample_dir"
        ((healed_flags+=1))
        ((skipped_complete+=1))
        continue
    fi

    if job_is_active "$sample_name"; then
        active_id="$(active_job_id "$sample_name" || true)"
        if [[ -n "$active_id" ]]; then
            log "[ACTIVE] $sample_name already has job $active_id; will wait before combining."
            dependency_job_ids+=("$active_id")
        else
            log "[WARNING] $sample_name appears active, but its JobID could not be resolved."
            ((failed_submission+=1))
        fi
        ((skipped_active+=1))
        continue
    fi

    # Incomplete/stale outputs are invalidated only when no equivalent job is active.
    rm -f "$success_flag" "$GLOBAL_SUCCESS_FLAG" "$COMBINE_SUCCESS_FLAG" "$COMBINE_WARNING_FLAG"

    {
        cat <<HEADER
#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH --job-name=${job_name}
#SBATCH --error=${job_stderr}
#SBATCH --output=${job_stdout}
#SBATCH --partition=k2-medpri,k2-bioinf,k2-lowpri,k2-himem
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
SAMPLE_NAME=$(printf '%q' "$sample_name")
CONTIGS=$(printf '%q' "$contigs")
SAMPLE_DIR=$(printf '%q' "$sample_dir")
OUTDIR=$(printf '%q' "$outdir")
DBNAME=$(printf '%q' "$DBNAME")
KRAKEN_ENV=$(printf '%q' "$KRAKEN_ENV")
SUCCESS_FLAG=$(printf '%q' "$success_flag")
HEADER
        cat <<'CHILD_BODY'

module load apps/anaconda3/2024.06/bin
source activate "$KRAKEN_ENV"
mkdir -p "$OUTDIR"

on_exit() {
    local rc="${1:-1}"
    trap - EXIT
    if [[ "$rc" -ne 0 ]]; then
        echo "ERROR: Kraken2 processing failed for ${SAMPLE_NAME} (exit $rc)." >&2
        rm -f "$SUCCESS_FLAG" \
          "${OUTDIR}/${SAMPLE_NAME}_kraken2.txt.tmp" \
          "${OUTDIR}/${SAMPLE_NAME}_kraken2_report.txt.tmp" \
          "${OUTDIR}/${SAMPLE_NAME}_kraken2_report_mpa.txt.tmp"
    fi
    exit "$rc"
}
trap 'on_exit "$?"' EXIT

[[ -s "$CONTIGS" ]] || { echo "ERROR: Contig FASTA missing: $CONTIGS" >&2; exit 1; }
[[ -d "$DBNAME" ]] || { echo "ERROR: Kraken2 database missing: $DBNAME" >&2; exit 1; }
for cmd in kraken2 seqtk python; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: Required command not found: $cmd" >&2; exit 1; }
done

FINAL_KRAKEN="${OUTDIR}/${SAMPLE_NAME}_kraken2.txt"
FINAL_REPORT="${OUTDIR}/${SAMPLE_NAME}_kraken2_report.txt"
FINAL_MPA="${OUTDIR}/${SAMPLE_NAME}_kraken2_report_mpa.txt"
TMP_KRAKEN="${FINAL_KRAKEN}.tmp"
TMP_REPORT="${FINAL_REPORT}.tmp"
TMP_MPA="${FINAL_MPA}.tmp"

generate_mpa_report() {
    python - "$1" "$2" <<'PYMPA'
import re
import sys

report_file, output_file = sys.argv[1], sys.argv[2]
prefix = {"D":"d", "K":"k", "P":"p", "C":"c", "O":"o", "F":"f", "G":"g", "S":"s"}
stack = []
with open(report_file, encoding="utf-8") as src, open(output_file, "w", encoding="utf-8") as dst:
    for line in src:
        raw = line.rstrip("\n")
        fields = raw.split("\t")
        if len(fields) < 6:
            continue
        try:
            clade = int(fields[1])
        except ValueError:
            continue
        rank = fields[3].strip()
        taxon_field = fields[5]
        depth = (len(taxon_field) - len(taxon_field.lstrip(" "))) // 2
        taxon = taxon_field.strip()
        taxon = re.sub(r"\s+\(taxid\s+\d+\)\s*$", "", taxon)
        taxon = taxon.replace("|", "").replace(".", "").replace(" ", "_")
        stack = stack[:depth]
        stack.append((rank, taxon))
        if rank not in prefix or clade <= 0 or not taxon:
            continue
        path = [f"{prefix[r]}__{name}" for r, name in stack if r in prefix and name]
        if path:
            dst.write("|".join(path) + "\t" + str(clade) + "\n")
PYMPA
}

# If classification and standard report are already current, repair only MPA
# when it is missing/stale. This avoids rerunning Kraken2 unnecessarily.
if [[ -s "$FINAL_KRAKEN" && -s "$FINAL_REPORT" ]] && \
   [[ "$CONTIGS" -ot "$FINAL_KRAKEN" ]] && \
   [[ "$CONTIGS" -ot "$FINAL_REPORT" ]]; then
    if [[ ! -s "$FINAL_MPA" || ! "$FINAL_REPORT" -ot "$FINAL_MPA" ]]; then
        echo "Repairing missing/stale MPA report only for ${SAMPLE_NAME}."
        rm -f "$TMP_MPA"
        generate_mpa_report "$FINAL_REPORT" "$TMP_MPA"
        [[ -s "$TMP_MPA" ]] || { echo "ERROR: MPA report repair failed." >&2; exit 1; }
        mv "$TMP_MPA" "$FINAL_MPA"
    fi
    [[ -s "$FINAL_MPA" && "$CONTIGS" -ot "$FINAL_MPA" && "$FINAL_REPORT" -ot "$FINAL_MPA" ]] || {
        echo "ERROR: Existing outputs are not current after MPA check." >&2
        exit 1
    }
    printf '%s\n' "$(date --iso-8601=seconds)" > "$SUCCESS_FLAG"
    echo "Existing current Kraken2 outputs verified/healed for ${SAMPLE_NAME}."
    trap - EXIT
    exit 0
fi

rm -f "$SUCCESS_FLAG" "$TMP_KRAKEN" "$TMP_REPORT" "$TMP_MPA"
echo "Running Kraken2 on contigs >= 2500 bp for ${SAMPLE_NAME}."
seqtk seq -L 2500 "$CONTIGS" | kraken2 \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --db "$DBNAME" \
    --use-names \
    --output "$TMP_KRAKEN" \
    --report "$TMP_REPORT" \
    /dev/stdin

[[ -s "$TMP_KRAKEN" ]] || { echo "ERROR: Kraken2 classification output is empty." >&2; exit 1; }
[[ -s "$TMP_REPORT" ]] || { echo "ERROR: Kraken2 standard report is empty." >&2; exit 1; }
generate_mpa_report "$TMP_REPORT" "$TMP_MPA"
[[ -s "$TMP_MPA" ]] || { echo "ERROR: MPA report generation failed." >&2; exit 1; }

mv "$TMP_KRAKEN" "$FINAL_KRAKEN"
mv "$TMP_REPORT" "$FINAL_REPORT"
mv "$TMP_MPA" "$FINAL_MPA"
[[ -s "$FINAL_KRAKEN" && -s "$FINAL_REPORT" && -s "$FINAL_MPA" ]] || { echo "ERROR: Final output validation failed." >&2; exit 1; }
[[ "$CONTIGS" -ot "$FINAL_KRAKEN" && "$CONTIGS" -ot "$FINAL_REPORT" && "$CONTIGS" -ot "$FINAL_MPA" ]] || { echo "ERROR: Final outputs are older than contigs." >&2; exit 1; }
printf '%s\n' "$(date --iso-8601=seconds)" > "$SUCCESS_FLAG"
echo "Kraken2 completed successfully for ${SAMPLE_NAME}; flag: $SUCCESS_FLAG"
CHILD_BODY
    } > "$job_script"
    chmod 750 "$job_script"

    if ! bash -n "$job_script"; then
        log "ERROR: Generated child job failed bash -n: $job_script"
        rm -f "$success_flag"
        ((failed_submission+=1))
        continue
    fi

    log "Submitting Kraken2 job for $sample_name..."
    if sbatch_output="$(sbatch "$job_script" 2>&1)"; then
        log "$sbatch_output"
        child_job_id="$(sed -n 's/^Submitted batch job \([0-9][0-9]*\).*$/\1/p' <<< "$sbatch_output" | head -n1)"
        if [[ -n "$child_job_id" ]]; then
            dependency_job_ids+=("$child_job_id")
            printf '%s\n' "$child_job_id" > "$jid_file"
        else
            log "WARNING: Could not parse submitted JobID for $sample_name; checker will still assess outputs."
        fi
        ((submitted+=1))
    else
        submit_rc=$?
        log "WARNING: sbatch failed for $sample_name (exit $submit_rc): $sbatch_output"
        rm -f "$success_flag" "$jid_file"
        ((failed_submission+=1))
    fi
done < "$SAMPLE_MANIFEST"

# ----------------------------------------------------------------------------
# Generate final checker. It waits for all submitted/active child jobs using
# afterany, then combines successful samples and writes warnings for failures.
# ----------------------------------------------------------------------------
{
    cat <<CHECK_HEADER
#!/bin/bash
#SBATCH --job-name=Kraken2_Contigs_Combined_Check
#SBATCH --error=${COMPLETION_ERR}
#SBATCH --output=${COMPLETION_LOG}
#SBATCH --partition=k2-hipri
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail
WORK_DIR=$(printf '%q' "$WORK_DIR")
JOB_DIR=$(printf '%q' "$JOB_DIR")
DBNAME=$(printf '%q' "$DBNAME")
KRAKEN_ENV=$(printf '%q' "$KRAKEN_ENV")
SAMPLE_MANIFEST=$(printf '%q' "$SAMPLE_MANIFEST")
OUTPUT_DIR=$(printf '%q' "$OUTPUT_DIR")
COMBINED_DIR=$(printf '%q' "$COMBINED_DIR")
SUMMARY_DIR=$(printf '%q' "$SUMMARY_DIR")
ABUNDANCE_DIR=$(printf '%q' "$ABUNDANCE_DIR")
VIS_DIR=$(printf '%q' "$VIS_DIR")
KRONA_DIR=$(printf '%q' "$KRONA_DIR")
GLOBAL_SUCCESS_FLAG=$(printf '%q' "$GLOBAL_SUCCESS_FLAG")
COMBINE_SUCCESS_FLAG=$(printf '%q' "$COMBINE_SUCCESS_FLAG")
COMBINE_WARNING_FLAG=$(printf '%q' "$COMBINE_WARNING_FLAG")
CHECK_HEADER
    cat <<'CHECK_BODY'

module load apps/anaconda3/2024.06/bin
source activate "$KRAKEN_ENV"
mkdir -p "$OUTPUT_DIR" "$COMBINED_DIR" "$SUMMARY_DIR" "$ABUNDANCE_DIR" "$VIS_DIR" "$KRONA_DIR"

on_exit() {
    local rc="${1:-1}"
    trap - EXIT
    if [[ "$rc" -ne 0 ]]; then
        rm -f "$COMBINE_SUCCESS_FLAG" "$GLOBAL_SUCCESS_FLAG"
        echo "ERROR: Kraken2 contig aggregation failed (exit $rc)." >&2
    fi
    exit "$rc"
}
trap 'on_exit "$?"' EXIT

STATUS_TSV="${SUMMARY_DIR}/Ran_Kraken2_contig_sample_status.tsv"
WARNINGS_TXT="${SUMMARY_DIR}/Ran_Kraken2_contig_warnings.txt"
SUCCESS_MANIFEST="${SUMMARY_DIR}/Ran_Kraken2_successful_contig_samples.tsv"

rm -f "$COMBINE_SUCCESS_FLAG" "$COMBINE_WARNING_FLAG" "$GLOBAL_SUCCESS_FLAG"
printf 'Sample\tStatus\tSample_directory\tContig_file\tSuccess_flag\tReason\n' > "$STATUS_TSV"
printf 'Kraken2 contig aggregation warnings\nGenerated: %s\n\n' "$(date --iso-8601=seconds)" > "$WARNINGS_TXT"
printf 'Sample\tSample_directory\tContig_file\tSuccess_flag\n' > "$SUCCESS_MANIFEST"

success_count=0
failure_count=0
total_count=0

while IFS=$'\t' read -r sample sample_dir contigs flag; do
    [[ "$sample" == "Sample" || -z "$sample" ]] && continue
    ((total_count+=1))
    outdir="${sample_dir}/${sample}_kraken2"
    kraken_file="${outdir}/${sample}_kraken2.txt"
    report_file="${outdir}/${sample}_kraken2_report.txt"
    mpa_file="${outdir}/${sample}_kraken2_report_mpa.txt"
    reason=""

    if [[ ! -s "$contigs" ]]; then
        reason="Missing or empty contigs FASTA"
    elif [[ ! -s "$kraken_file" ]]; then
        reason="Missing/empty Kraken2 classification output"
    elif [[ ! -s "$report_file" ]]; then
        reason="Missing/empty Kraken2 standard report"
    elif [[ ! -s "$mpa_file" ]]; then
        reason="Missing/empty MPA report"
    elif [[ ! "$contigs" -ot "$kraken_file" || ! "$contigs" -ot "$report_file" || ! "$contigs" -ot "$mpa_file" ]]; then
        reason="One or more outputs are older than the contig FASTA"
    elif [[ ! "$report_file" -ot "$mpa_file" ]]; then
        reason="MPA report is older than the standard report"
    fi

    if [[ -z "$reason" ]]; then
        if [[ ! -f "$flag" ]]; then
            printf '%s\n' "$(date --iso-8601=seconds)" > "$flag"
            echo "[HEAL] Recreated success flag for $sample: $flag"
        fi
        printf '%s\tSUCCESS\t%s\t%s\t%s\t\n' "$sample" "$sample_dir" "$contigs" "$flag" >> "$STATUS_TSV"
        printf '%s\t%s\t%s\t%s\n' "$sample" "$sample_dir" "$contigs" "$flag" >> "$SUCCESS_MANIFEST"
        ((success_count+=1))
    else
        rm -f "$flag"
        printf '%s\tFAILED_OR_INCOMPLETE\t%s\t%s\t%s\t%s\n' "$sample" "$sample_dir" "$contigs" "$flag" "$reason" >> "$STATUS_TSV"
        printf 'WARNING: %s - %s\n' "$sample" "$reason" >> "$WARNINGS_TXT"
        echo "[WARNING] $sample: $reason" >&2
        ((failure_count+=1))
    fi
done < "$SAMPLE_MANIFEST"

if [[ "$total_count" -eq 0 ]]; then
    echo "ERROR: No samples were found in the sample manifest." >&2
    exit 1
fi

printf 'Samples checked: %s\nSuccessful samples: %s\nFailed/incomplete samples: %s\n\n' \
  "$total_count" "$success_count" "$failure_count" >> "$WARNINGS_TXT"

# All aggregation logic is embedded here: no separately maintained .py file.
python - "$SUCCESS_MANIFEST" "$STATUS_TSV" "$OUTPUT_DIR" "$WARNINGS_TXT" <<'PYAGG'
import csv
import gzip
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict

success_manifest, status_file, out_root, warnings_file = sys.argv[1:]
combined_dir = os.path.join(out_root, "combined")
summary_dir = os.path.join(out_root, "summaries")
abundance_root = os.path.join(out_root, "abundance")
vis_root = os.path.join(out_root, "visualisations")
krona_root = os.path.join(out_root, "Krona")
for path in (combined_dir, summary_dir, abundance_root, vis_root, krona_root):
    os.makedirs(path, exist_ok=True)

with open(success_manifest, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    successes = list(reader)
with open(status_file, newline="", encoding="utf-8") as f:
    status_rows = list(csv.DictReader(f, delimiter="\t"))

samples = [row["Sample"] for row in successes]
status_by_sample = {row["Sample"]: row for row in status_rows}
rank_map = {"D":"domain", "K":"kingdom", "P":"phylum", "C":"class", "O":"order", "F":"family", "G":"genus", "S":"species"}

out_contigs = os.path.join(combined_dir, "Ran_Kraken2_combined_contig_level.tsv.gz")
out_reports = os.path.join(combined_dir, "Ran_Kraken2_combined_report.tsv.gz")
out_mpa = os.path.join(combined_dir, "Ran_Kraken2_combined_mpa.tsv.gz")

# Combine sequence-level classification records, retaining explicit sample IDs.
summary_values = {}
with gzip.open(out_contigs, "wt", encoding="utf-8", newline="") as out:
    out.write("Sample\tKraken_classification\tTaxID\tContig_ID\tLength\tLCA_mapping\n")
    for row in successes:
        sample = row["Sample"]
        outdir = os.path.join(row["Sample_directory"], f"{sample}_kraken2")
        path = os.path.join(outdir, f"{sample}_kraken2.txt")
        total = classified = unclassified = malformed = 0
        with open(path, encoding="utf-8", errors="replace") as src:
            for line in src:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    malformed += 1
                    continue
                cls, contig_id, taxid, length, lca = fields[:5]
                out.write("\t".join((sample, cls, taxid, contig_id, length, lca)) + "\n")
                total += 1
                if cls == "C": classified += 1
                elif cls == "U": unclassified += 1
        summary_values[sample] = {"total": total, "classified": classified, "unclassified": unclassified}
        if malformed:
            with open(warnings_file, "a", encoding="utf-8") as warn:
                warn.write(f"WARNING: {sample} classification output had {malformed} malformed rows.\n")

# Combine standard report rows and MPA rows for successful samples only.
with gzip.open(out_reports, "wt", encoding="utf-8", newline="") as out, \
     gzip.open(out_mpa, "wt", encoding="utf-8", newline="") as mpa_out:
    out.write("Sample\tPercentage\tClade_contigs\tTaxon_contigs\tRank\tTaxID\tName\n")
    mpa_out.write("Sample\tTaxonomy_path\tClade_contigs\n")
    for row in successes:
        sample = row["Sample"]
        outdir = os.path.join(row["Sample_directory"], f"{sample}_kraken2")
        report_path = os.path.join(outdir, f"{sample}_kraken2_report.txt")
        with open(report_path, encoding="utf-8", errors="replace") as src:
            for line in src:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                name = fields[5].lstrip(" \t")
                out.write("\t".join((sample, fields[0], fields[1], fields[2], fields[3], fields[4], name)) + "\n")
        mpa_path = os.path.join(outdir, f"{sample}_kraken2_report_mpa.txt")
        with open(mpa_path, encoding="utf-8", errors="replace") as src:
            for line in src:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    mpa_out.write(f"{sample}\t{fields[0]}\t{fields[1]}\n")

# Sample-level summary includes failures explicitly; failed samples get NA, not zero.
summary_file = os.path.join(summary_dir, "Ran_Kraken2_contig_sample_summary.tsv")
with open(summary_file, "w", encoding="utf-8", newline="") as out:
    out.write("Sample\tStatus\tTotal_contigs_ge_2500bp\tClassified_contigs\tUnclassified_contigs\tPercent_classified\tPercent_unclassified\tNote\n")
    for status in status_rows:
        sample = status["Sample"]
        if status["Status"] == "SUCCESS" and sample in summary_values:
            vals = summary_values[sample]
            total = vals["total"]
            classified = vals["classified"]
            unclassified = vals["unclassified"]
            pc = 100.0 * classified / total if total else 0.0
            pu = 100.0 * unclassified / total if total else 0.0
            note = ""
            out.write(f"{sample}\tSUCCESS\t{total}\t{classified}\t{unclassified}\t{pc:.4f}\t{pu:.4f}\t{note}\n")
        else:
            out.write(f"{sample}\tFAILED_OR_INCOMPLETE\tNA\tNA\tNA\tNA\tNA\t{status.get('Reason','')}\n")

# Abundance tables use only successful samples; failures are visible in the status table.
counts = {rank: defaultdict(dict) for rank in rank_map.values()}
names = {rank: {} for rank in rank_map.values()}
for row in successes:
    sample = row["Sample"]
    outdir = os.path.join(row["Sample_directory"], f"{sample}_kraken2")
    report_path = os.path.join(outdir, f"{sample}_kraken2_report.txt")
    with open(report_path, encoding="utf-8", errors="replace") as src:
        for line in src:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[3] not in rank_map:
                continue
            rank = rank_map[fields[3]]
            try: clade = int(fields[1])
            except ValueError: continue
            taxid = fields[4]
            name = fields[5].strip()
            counts[rank].setdefault(taxid, {})[sample] = clade
            names[rank][taxid] = name

sample_totals = {s: summary_values.get(s, {}).get("total", 0) for s in samples}
for rank in ("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"):
    rank_dir = os.path.join(abundance_root, rank)
    os.makedirs(rank_dir, exist_ok=True)
    rank_counts = counts[rank]
    totals = {taxid: sum(values.values()) for taxid, values in rank_counts.items()}
    ordered = sorted(totals, key=lambda t: (-totals[t], names[rank].get(t, ""), t))
    count_path = os.path.join(rank_dir, f"Ran_Kraken2_contigs_{rank}_counts.tsv")
    abundance_path = os.path.join(rank_dir, f"Ran_Kraken2_contigs_{rank}_relative_abundance.tsv")
    presence_path = os.path.join(rank_dir, f"Ran_Kraken2_contigs_{rank}_presence.tsv")
    metadata_path = os.path.join(rank_dir, f"Ran_Kraken2_contigs_{rank}_taxa_metadata.tsv")
    with open(count_path, "w", encoding="utf-8") as out:
        out.write("TaxID\tTaxon\t" + "\t".join(samples) + "\tTotal\n")
        for taxid in ordered:
            vals = [rank_counts[taxid].get(s, 0) for s in samples]
            out.write(f"{taxid}\t{names[rank].get(taxid,'')}\t" + "\t".join(map(str, vals)) + f"\t{sum(vals)}\n")
    with open(abundance_path, "w", encoding="utf-8") as out:
        out.write("TaxID\tTaxon\t" + "\t".join(samples) + "\tMean_relative_abundance\n")
        for taxid in ordered:
            vals = [(100.0 * rank_counts[taxid].get(s, 0) / sample_totals[s]) if sample_totals[s] else 0.0 for s in samples]
            mean = sum(vals) / len(vals) if vals else 0.0
            out.write(f"{taxid}\t{names[rank].get(taxid,'')}\t" + "\t".join(f"{x:.6f}" for x in vals) + f"\t{mean:.6f}\n")
    with open(presence_path, "w", encoding="utf-8") as out:
        out.write("TaxID\tTaxon\t" + "\t".join(samples) + "\n")
        for taxid in ordered:
            vals = [1 if rank_counts[taxid].get(s, 0) > 0 else 0 for s in samples]
            out.write(f"{taxid}\t{names[rank].get(taxid,'')}\t" + "\t".join(map(str, vals)) + "\n")
    with open(metadata_path, "w", encoding="utf-8") as out:
        out.write("TaxID\tTaxon\tRank\tTotal_contigs_across_samples\n")
        for taxid in ordered:
            out.write(f"{taxid}\t{names[rank].get(taxid,'')}\t{rank}\t{totals[taxid]}\n")

# Generate plots when matplotlib is available. Failures here are warnings and do
# not discard the combined classification/taxonomy tables.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    successful_summary = [s for s in samples if s in summary_values]
    if successful_summary:
        rows = [summary_values[s] for s in successful_summary]
        classified_pct = [100.0 * r["classified"] / r["total"] if r["total"] else 0.0 for r in rows]
        unclassified_pct = [100.0 * r["unclassified"] / r["total"] if r["total"] else 0.0 for r in rows]
        totals = [r["total"] for r in rows]
        fig, ax = plt.subplots(figsize=(max(12, len(successful_summary)*0.22), 7))
        ax.bar(range(len(successful_summary)), classified_pct, label="Classified")
        ax.bar(range(len(successful_summary)), unclassified_pct, bottom=classified_pct, label="Unclassified")
        ax.set_ylabel("Percentage of contigs >= 2500 bp")
        ax.set_xticks(range(len(successful_summary))); ax.set_xticklabels(successful_summary, rotation=90, fontsize=7)
        ax.set_ylim(0, 100); ax.set_title("Kraken2 contig classification rate"); ax.legend(); fig.tight_layout()
        fig.savefig(os.path.join(vis_root, "classified_contigs_per_sample.png"), dpi=180, bbox_inches="tight"); plt.close(fig)
        order = sorted(range(len(totals)), key=lambda i: totals[i])
        fig, ax = plt.subplots(figsize=(max(12, len(successful_summary)*0.22), 7))
        ax.bar(range(len(order)), [totals[i] for i in order])
        ax.set_ylabel("Contigs classified by Kraken2 (>=2500 bp)")
        ax.set_xticks(range(len(order))); ax.set_xticklabels([successful_summary[i] for i in order], rotation=90, fontsize=7)
        ax.set_title("Number of contigs per sample"); fig.tight_layout()
        fig.savefig(os.path.join(vis_root, "total_contigs_per_sample.png"), dpi=180, bbox_inches="tight"); plt.close(fig)

        def stacked(rank, top_n, filename, title):
            path = os.path.join(abundance_root, rank, f"Ran_Kraken2_contigs_{rank}_relative_abundance.tsv")
            with open(path, newline="", encoding="utf-8") as f:
                recs = list(csv.DictReader(f, delimiter="\t"))
            ranked = sorted(recs, key=lambda r: float(r["Mean_relative_abundance"]), reverse=True)[:top_n]
            xs = list(range(len(successful_summary))); bottom = [0.0] * len(xs)
            fig, ax = plt.subplots(figsize=(max(14, len(xs)*0.24), 8))
            top_ids = {r["TaxID"] for r in ranked}
            for rec in ranked:
                vals = [float(rec.get(s, 0) or 0) for s in successful_summary]
                ax.bar(xs, vals, bottom=bottom, label=rec["Taxon"])
                bottom = [a+b for a,b in zip(bottom, vals)]
            other = []
            for s in successful_summary:
                other.append(sum(float(r.get(s,0) or 0) for r in recs if r["TaxID"] not in top_ids))
            ax.bar(xs, other, bottom=bottom, label="Other")
            ax.set_ylabel("Relative abundance (%)"); ax.set_xlabel("Sample"); ax.set_title(title)
            ax.set_xticks(xs); ax.set_xticklabels(successful_summary, rotation=90, fontsize=7)
            if ranked: ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
            fig.tight_layout(); fig.savefig(os.path.join(vis_root, filename), dpi=180, bbox_inches="tight"); plt.close(fig)
        stacked("genus", 20, "top_20_genera_stacked_barplot.png", "Top 20 genera across successful samples")
        stacked("family", 20, "top_20_families_stacked_barplot.png", "Top 20 families across successful samples")
        stacked("phylum", 15, "top_15_phyla_stacked_barplot.png", "Top 15 phyla across successful samples")

        genus_path = os.path.join(abundance_root, "genus", "Ran_Kraken2_contigs_genus_relative_abundance.tsv")
        with open(genus_path, newline="", encoding="utf-8") as f: genus_rows = list(csv.DictReader(f, delimiter="\t"))
        top_genus = sorted(genus_rows, key=lambda r: float(r["Mean_relative_abundance"]), reverse=True)[:30]
        if top_genus:
            matrix = [[float(r.get(s,0) or 0) for s in successful_summary] for r in top_genus]
            fig, ax = plt.subplots(figsize=(max(12, len(successful_summary)*0.25), max(7, len(top_genus)*0.28)))
            im = ax.imshow(matrix, aspect="auto", interpolation="nearest")
            ax.set_yticks(range(len(top_genus))); ax.set_yticklabels([r["Taxon"] for r in top_genus], fontsize=7)
            ax.set_xticks(range(len(successful_summary))); ax.set_xticklabels(successful_summary, rotation=90, fontsize=6)
            ax.set_title("Top 30 genera by mean relative abundance"); fig.colorbar(im, ax=ax, label="Relative abundance (%)"); fig.tight_layout()
            fig.savefig(os.path.join(vis_root, "top_30_genera_heatmap.png"), dpi=180, bbox_inches="tight"); plt.close(fig)
        richness = [sum(1 for r in genus_rows if float(r.get(s,0) or 0) > 0) for s in successful_summary]
        order = sorted(range(len(richness)), key=lambda i: richness[i])
        fig, ax = plt.subplots(figsize=(max(12, len(successful_summary)*0.22), 7))
        ax.bar(range(len(order)), [richness[i] for i in order])
        ax.set_ylabel("Number of detected genera"); ax.set_xticks(range(len(order)))
        ax.set_xticklabels([successful_summary[i] for i in order], rotation=90, fontsize=7)
        ax.set_title("Kraken2 genus-level richness"); fig.tight_layout()
        fig.savefig(os.path.join(vis_root, "genus_richness_per_sample.png"), dpi=180, bbox_inches="tight"); plt.close(fig)
    else:
        with open(os.path.join(vis_root, "NO_SUCCESSFUL_SAMPLES.txt"), "w", encoding="utf-8") as f:
            f.write("No successful Kraken2 contig samples were available for plotting. See summaries/Ran_Kraken2_contig_sample_status.tsv.\n")
except Exception as exc:
    with open(warnings_file, "a", encoding="utf-8") as warn:
        warn.write(f"WARNING: Visualisation generation failed or was unavailable: {exc}\n")

# Krona is optional. Use the per-sample embedded MPA-style taxonomy paths.
kt = shutil.which("ktImportText")
if kt:
    for row in successes:
        sample = row["Sample"]
        outdir = os.path.join(row["Sample_directory"], f"{sample}_kraken2")
        mpa_path = os.path.join(outdir, f"{sample}_kraken2_report_mpa.txt")
        input_path = os.path.join(krona_root, f"{sample}_Krona_input.tsv")
        html_path = os.path.join(krona_root, f"{sample}_Krona.html")
        with open(mpa_path, encoding="utf-8", errors="replace") as src, open(input_path, "w", encoding="utf-8") as dst:
            for line in src:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2: continue
                try: n = int(fields[1])
                except ValueError: continue
                names = [p.split("__",1)[-1] for p in fields[0].split("|")]
                dst.write(str(n) + "\t" + "\t".join(names) + "\n")
        try:
            subprocess.run([kt, "-o", html_path, input_path], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        except Exception as exc:
            with open(warnings_file, "a", encoding="utf-8") as warn:
                warn.write(f"WARNING: Krona generation failed for {sample}: {exc}\n")
else:
    with open(warnings_file, "a", encoding="utf-8") as warn:
        warn.write("NOTE: ktImportText not found; Krona HTML generation was skipped.\n")

print(f"Aggregated {len(successes)} successful samples; {len(status_rows)-len(successes)} sample(s) failed or were incomplete.")
PYAGG

# Validate core aggregate outputs before publishing a completion state.
for f in \
  "${COMBINED_DIR}/Ran_Kraken2_combined_contig_level.tsv.gz" \
  "${COMBINED_DIR}/Ran_Kraken2_combined_report.tsv.gz" \
  "${COMBINED_DIR}/Ran_Kraken2_combined_mpa.tsv.gz" \
  "${SUMMARY_DIR}/Ran_Kraken2_contig_sample_status.tsv" \
  "${SUMMARY_DIR}/Ran_Kraken2_contig_sample_summary.tsv" \
  "${SUMMARY_DIR}/Ran_Kraken2_contig_warnings.txt"; do
    [[ -s "$f" ]] || { echo "ERROR: Required combined output missing/empty: $f" >&2; exit 1; }
done

if [[ "$failure_count" -eq 0 && "$success_count" -gt 0 ]]; then
    printf '%s\n' "$(date --iso-8601=seconds)" > "$COMBINE_SUCCESS_FLAG"
    printf '%s\n' "$(date --iso-8601=seconds)" > "$GLOBAL_SUCCESS_FLAG"
    rm -f "$COMBINE_WARNING_FLAG"
    echo "Kraken2 contig aggregation completed without sample failures."
    echo "Success flags: $COMBINE_SUCCESS_FLAG ; $GLOBAL_SUCCESS_FLAG"
else
    printf '%s\n' "$(date --iso-8601=seconds)" > "$COMBINE_WARNING_FLAG"
    rm -f "$COMBINE_SUCCESS_FLAG" "$GLOBAL_SUCCESS_FLAG"
    echo "Kraken2 contig aggregation completed with warnings." >&2
    echo "Successful samples: $success_count / $total_count" >&2
    echo "Details: $WARNINGS_TXT" >&2
fi

echo "Combined output directory: $OUTPUT_DIR"
trap - EXIT
exit 0
CHECK_BODY
} > "$COMPLETION_JOB"
chmod 750 "$COMPLETION_JOB"
if ! bash -n "$COMPLETION_JOB"; then
    log "ERROR: Generated completion checker failed bash -n: $COMPLETION_JOB"
    exit 1
fi

# Include just-submitted and already-active jobs; afterany is intentional so
# one failed sample cannot prevent aggregation of successful samples.
if ((${#dependency_job_ids[@]} > 0)); then
    declare -A seen_ids=()
    unique_ids=()
    for jid in "${dependency_job_ids[@]}"; do
        if [[ -n "$jid" && -z "${seen_ids[$jid]+x}" ]]; then
            seen_ids[$jid]=1
            unique_ids+=("$jid")
        fi
    done
    dependency_list="$(IFS=:; echo "${unique_ids[*]}")"
    log "Submitting combined-results checker with afterany dependency: $dependency_list"
    if checker_output="$(sbatch --dependency="afterany:${dependency_list}" "$COMPLETION_JOB" 2>&1)"; then
        log "$checker_output"
    else
        checker_rc=$?
        log "ERROR: Could not submit combined-results checker (exit $checker_rc): $checker_output"
        exit 1
    fi
else
    log "No active/new Kraken2 jobs to wait for; submitting checker to combine existing sample outputs."
    if checker_output="$(sbatch "$COMPLETION_JOB" 2>&1)"; then
        log "$checker_output"
    else
        checker_rc=$?
        log "ERROR: Could not submit combined-results checker (exit $checker_rc): $checker_output"
        exit 1
    fi
fi

printf '%s\n' "$(date --iso-8601=seconds)" > "$PREP_SUCCESS_FLAG"
log "Preparation/submission flag: $PREP_SUCCESS_FLAG"
log "Sample-level flags: <sample>/<sample>_kraken2.success"
log "Combined output directory: $OUTPUT_DIR"
log "Preparation complete; combination runs after all active/submitted jobs finish."
trap - EXIT
exit 0