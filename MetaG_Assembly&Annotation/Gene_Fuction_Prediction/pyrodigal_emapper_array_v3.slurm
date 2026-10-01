#!/bin/bash

#SBATCH --job-name=Prepare_Pyro_Eggnog_Ran
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --error=/users/3057556/jobs/Prepare_Pyro_Eggnog_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_Pyro_Eggnog_Ran_%j.log
#SBATCH --partition=k2-medpri,k2-bioinf,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# =========================================================================
# Combined Pyrodigal + eggNOG-mapper preparation/submission script
#
# Workflow:
#
#   MetaSPAdes assembly
#          |
#          | assembly QC
#          v
#      Pyrodigal
#          |
#          | protein/output QC
#          v
#      eggNOG-mapper
#
# One Pyrodigal SLURM job and one eggNOG-mapper SLURM job are generated
# per sample.
#
# IMPORTANT:
#   - This is a preparation/submission script.
#   - The actual Pyrodigal and eggNOG jobs run independently.
#   - eggNOG is submitted with --dependency=afterok:<Pyrodigal_JOBID>
#     so that it cannot start until Pyrodigal has successfully completed.
#
# Restartable:
#   - Completed Pyrodigal samples are skipped.
#   - Completed eggNOG samples are skipped.
#   - Existing incomplete outputs are removed only for the stage being
#     rerun.
#   - Assembly QC is repeated because it is cheap and protects against
#     accidentally processing a changed/corrupt assembly.
#
# =========================================================================

set -euo pipefail

# =========================================================================
# 1. Load Anaconda
# =========================================================================

module load apps/anaconda3/2024.06/bin

# =========================================================================
# 2. Define paths
# =========================================================================

WORK_DIR="/mnt/scratch2/users/3057556/Ran/40-1310848928/00_fastq"

JOB_DIR="/users/3057556/jobs/Ran_pyro_eggnog_jobs"

PYRODIGAL_ENV="/mnt/scratch2/igfs-anaconda/conda-envs/pyrodigal_3.6.3"

EGGNOG_ENV="/mnt/scratch2/igfs-anaconda/conda-envs/eggnog-mapper-2.1.12"

EGGNOG_DB="/mnt/scratch2/igfs-anaconda/conda-envs/eggnog-mapper-2.1.12/db"

# -------------------------------------------------------------------------
# Directory tags
# -------------------------------------------------------------------------

DIR_TAGS=("4" "N" "P")

# -------------------------------------------------------------------------
# Resource requests
# -------------------------------------------------------------------------

PYRO_CPUS=10
PYRO_MEM=100G
PYRO_TIME="03:00:00"

EGGNOG_CPUS=20
EGGNOG_MEM=200G
EGGNOG_TIME="24:00:00"

# -------------------------------------------------------------------------
# Minimum assembly size warning threshold
#
# This is deliberately a WARNING rather than an automatic failure.
# Metagenomic assemblies can legitimately vary substantially.
# -------------------------------------------------------------------------

MIN_ASSEMBLY_BP_WARNING=10000

# =========================================================================
# 3. Create job directory
# =========================================================================

mkdir -p "$JOB_DIR"

cd "$WORK_DIR" || {
    echo "ERROR: Could not change to:"
    echo "$WORK_DIR"
    exit 1
}

echo
echo "============================================================"
echo "Combined Pyrodigal + eggNOG-mapper preparation"
echo "============================================================"
echo
echo "Working directory:"
echo "  $WORK_DIR"
echo
echo "Job directory:"
echo "  $JOB_DIR"
echo

# =========================================================================
# 4. Function: assembly QC
#
# This is run by the preparation script before a Pyrodigal job is made.
#
# It creates:
#
#   ${sample}_assembly_qc.tsv
#
# This contains basic assembly statistics.
# =========================================================================

check_assembly() {

    local contigs="$1"
    local sample="$2"
    local qc_file="$3"

    echo
    echo "------------------------------------------------------------"
    echo "Assembly QC: $sample"
    echo "------------------------------------------------------------"

    if [[ ! -f "$contigs" ]]; then

        echo "ERROR: Assembly does not exist:"
        echo "  $contigs"

        return 1

    fi

    if [[ ! -s "$contigs" ]]; then

        echo "ERROR: Assembly is empty:"
        echo "  $contigs"

        return 1

    fi

    # ---------------------------------------------------------------------
    # Parse FASTA and perform basic validation.
    #
    # Python is used here because it is considerably safer than trying
    # to perform complete FASTA validation using awk/sed.
    # ---------------------------------------------------------------------

    python - "$contigs" "$qc_file" "$MIN_ASSEMBLY_BP_WARNING" <<'PY'

import sys
from collections import Counter

fasta = sys.argv[1]
qc_file = sys.argv[2]
min_warning = int(sys.argv[3])

valid_bases = set("ACGTNacgtn")

sequence_count = 0
total_bp = 0
invalid_chars = Counter()
duplicate_ids = []
seen_ids = set()
empty_sequences = 0
lengths = []

current_id = None
current_seq = []

def process_sequence(seq_id, seq):

    global sequence_count
    global total_bp
    global empty_sequences

    if seq_id is None:
        return

    if not seq:
        empty_sequences += 1
        return

    sequence_count += 1

    length = len(seq)

    lengths.append(length)
    total_bp += length

    if seq_id in seen_ids:
        duplicate_ids.append(seq_id)

    seen_ids.add(seq_id)

    for base in seq:
        if base not in valid_bases:
            invalid_chars[base] += 1


with open(fasta, "r") as f:

    for line_number, raw in enumerate(f, 1):

        line = raw.rstrip("\r\n")

        if not line:
            continue

        if line.startswith(">"):

            process_sequence(
                current_id,
                "".join(current_seq)
            )

            header = line[1:].strip()

            if not header:
                print(
                    f"ERROR: Empty FASTA header at line {line_number}",
                    file=sys.stderr
                )
                sys.exit(1)

            # FASTA identifier = first whitespace-delimited token
            current_id = header.split()[0]

            current_seq = []

        else:

            if current_id is None:

                print(
                    f"ERROR: Sequence data encountered before "
                    f"first FASTA header at line {line_number}",
                    file=sys.stderr
                )

                sys.exit(1)

            current_seq.append(
                "".join(line.split())
            )

process_sequence(
    current_id,
    "".join(current_seq)
)

if sequence_count == 0:

    print(
        "ERROR: No sequences found in FASTA.",
        file=sys.stderr
    )

    sys.exit(1)

if empty_sequences > 0:

    print(
        f"ERROR: {empty_sequences} empty FASTA sequences found.",
        file=sys.stderr
    )

    sys.exit(1)

if invalid_chars:

    print(
        "ERROR: Invalid characters found in assembly:",
        file=sys.stderr
    )

    for char, count in sorted(invalid_chars.items()):

        print(
            f"  {repr(char)}: {count}",
            file=sys.stderr
        )

    sys.exit(1)

# -------------------------------------------------------------------------
# Assembly statistics
# -------------------------------------------------------------------------

lengths.sort(reverse=True)

largest = lengths[0]

# N50
half = total_bp / 2
running = 0
n50 = 0

for length in lengths:

    running += length

    if running >= half:

        n50 = length
        break

# N90
target = total_bp * 0.90
running = 0
n90 = 0

for length in lengths:

    running += length

    if running >= target:

        n90 = length
        break

mean_length = total_bp / sequence_count

gc_bases = 0

# Re-read sequence to calculate GC
with open(fasta, "r") as f:

    sequence = []

    for raw in f:

        line = raw.strip()

        if not line:
            continue

        if line.startswith(">"):

            if sequence:

                seq = "".join(sequence)

                gc_bases += sum(
                    1
                    for base in seq
                    if base in "GCgc"
                )

            sequence = []

        else:

            sequence.append(
                "".join(line.split())
            )

    if sequence:

        seq = "".join(sequence)

        gc_bases += sum(
            1
            for base in seq
            if base in "GCgc"
        )

gc_percent = (
    (gc_bases / total_bp) * 100
    if total_bp > 0
    else 0
)

warning = ""

if total_bp < min_warning:

    warning = (
        f"WARNING: assembly is only {total_bp} bp"
    )

# -------------------------------------------------------------------------
# Write QC file
# -------------------------------------------------------------------------

with open(qc_file, "w") as out:

    out.write(
        "Sample\tContigs\tTotal_bp\tLargest_contig_bp\t"
        "N50_bp\tN90_bp\tMean_contig_bp\tGC_percent\t"
        "Duplicate_IDs\tWarning\n"
    )

    out.write(
        f"{sample}\t"
        f"{sequence_count}\t"
        f"{total_bp}\t"
        f"{largest}\t"
        f"{n50}\t"
        f"{n90}\t"
        f"{mean_length:.2f}\t"
        f"{gc_percent:.4f}\t"
        f"{len(duplicate_ids)}\t"
        f"{warning}\n"
    )

print()
print("Assembly QC:")
print(f"  Contigs:           {sequence_count:,}")
print(f"  Total bp:          {total_bp:,}")
print(f"  Largest contig:    {largest:,} bp")
print(f"  N50:               {n50:,} bp")
print(f"  N90:               {n90:,} bp")
print(f"  Mean contig:       {mean_length:,.2f} bp")
print(f"  GC:                {gc_percent:.4f}%")
print(f"  Duplicate IDs:     {len(duplicate_ids)}")

if warning:
    print()
    print(warning)

PY

    if [[ ! -s "$qc_file" ]]; then

        echo "ERROR: Assembly QC output was not created."

        return 1

    fi

    echo "Assembly QC passed."

    return 0
}


# =========================================================================
# 5. Loop through samples
# =========================================================================

for TAG in "${DIR_TAGS[@]}"; do

    echo
    echo "============================================================"
    echo "Searching directories beginning with: $TAG"
    echo "============================================================"

    for dir in "$TAG"*; do

        [[ -d "$dir" ]] || continue

        sample_name=$(basename "$dir")

        echo
        echo "############################################################"
        echo "Sample: $sample_name"
        echo "############################################################"

        # =================================================================
        # 5A. Locate MetaSPAdes assembly
        # =================================================================

        contigs_file=$(find "$dir" \
            -maxdepth 3 \
            -type f \
            -name "contigs.fasta" \
            -path "*_metaspades/*" \
            -print -quit)

        if [[ -z "$contigs_file" ]]; then

            echo "WARNING: No MetaSPAdes contigs.fasta found."
            echo "Skipping $sample_name."

            continue

        fi

        echo "Assembly:"
        echo "  $contigs_file"

        # =================================================================
        # 5B. Require MetaSPAdes completion flag
        # =================================================================

        metaspades_flag="${dir}/${sample_name}_metaspades_assembly_done.flag"

        if [[ ! -f "$metaspades_flag" ]]; then

            echo
            echo "WARNING: MetaSPAdes completion flag missing:"
            echo "  $metaspades_flag"
            echo
            echo "The assembly exists but MetaSPAdes did not leave"
            echo "its completion flag. Skipping."

            continue

        fi

        # =================================================================
        # 5C. Assembly QC
        # =================================================================

        ASSEMBLY_QC="${dir}/${sample_name}_assembly_qc.tsv"

        if ! check_assembly \
            "$contigs_file" \
            "$sample_name" \
            "$ASSEMBLY_QC"
        then

            echo
            echo "ERROR: Assembly QC failed."
            echo "Pyrodigal will NOT be submitted for $sample_name."

            continue

        fi

        # =================================================================
        # 5D. Define Pyrodigal outputs
        # =================================================================

        PYRO_DIR="${dir}/${sample_name}_pyrodigal"

        PYRO_GFF="${PYRO_DIR}/${sample_name}_pyrodigal.gff"
        PYRO_NT="${PYRO_DIR}/${sample_name}_pyrodigal_nt.fa"
        PYRO_AA="${PYRO_DIR}/${sample_name}_pyrodigal_aa.fa"

        PYRO_FLAG="${PYRO_DIR}/${sample_name}_pyrodigal_done.flag"

        PYRO_QC="${PYRO_DIR}/${sample_name}_pyrodigal_qc.tsv"

        PYRO_JOB="${JOB_DIR}/${sample_name}_pyrodigal_job.sh"

        # =================================================================
        # 5E. Determine whether Pyrodigal is already complete
        # =================================================================

        if [[ -f "$PYRO_FLAG" \
              && -s "$PYRO_GFF" \
              && -s "$PYRO_NT" \
              && -s "$PYRO_AA" \
              && -s "$PYRO_QC" ]]; then

            echo
            echo "Pyrodigal already completed:"
            echo "  $PYRO_FLAG"

            PYRO_COMPLETE=1

        elif [[ -s "$PYRO_GFF" \
                && -s "$PYRO_NT" \
                && -s "$PYRO_AA" \
                && -s "$PYRO_QC" ]]; then

            echo
            echo "Pyrodigal outputs exist but completion flag is missing."

            echo "Treating outputs as complete because all validated"
            echo "outputs are present."

            printf "%s\n" \
                "$(date --iso-8601=seconds)" \
                > "$PYRO_FLAG"

            PYRO_COMPLETE=1

        else

            PYRO_COMPLETE=0

        fi

        # =================================================================
        # 5F. Define eggNOG outputs
        # =================================================================

        EGGNOG_DIR="${dir}/${sample_name}_eggnog_mapper"

        EGGNOG_PREFIX="${sample_name}_pyrodigal_eggnog_mapped"

        EGGNOG_XLSX="${EGGNOG_DIR}/${EGGNOG_PREFIX}.emapper.annotations.xlsx"

        EGGNOG_ANNOTATIONS="${EGGNOG_DIR}/${EGGNOG_PREFIX}.emapper.annotations"

        EGGNOG_FLAG="${EGGNOG_DIR}/${sample_name}_eggnogmapper.success"

        EGGNOG_QC="${EGGNOG_DIR}/${sample_name}_eggnogmapper_qc.tsv"

        EGGNOG_JOB="${JOB_DIR}/${sample_name}_eggnogMapper_job.sh"

        # =================================================================
        # 5G. Check whether eggNOG is already complete
        # =================================================================

        if [[ -f "$EGGNOG_FLAG" \
              && -s "$EGGNOG_XLSX" \
              && -s "$EGGNOG_QC" ]]; then

            echo
            echo "eggNOG-mapper already completed:"
            echo "  $EGGNOG_FLAG"

            EGGNOG_COMPLETE=1

        else

            EGGNOG_COMPLETE=0

        fi

        # =================================================================
        # 5H. If both stages are complete, nothing to submit
        # =================================================================

        if [[ "$PYRO_COMPLETE" -eq 1 \
              && "$EGGNOG_COMPLETE" -eq 1 ]]; then

            echo
            echo "============================================================"
            echo "$sample_name is already fully processed."
            echo "============================================================"

            continue

        fi

        # =================================================================
        # 5I. Generate Pyrodigal job if necessary
        # =================================================================

        if [[ "$PYRO_COMPLETE" -eq 0 ]]; then

            echo
            echo "Preparing Pyrodigal job."

            # -------------------------------------------------------------
            # Remove incomplete Pyrodigal outputs.
            # -------------------------------------------------------------

            mkdir -p "$PYRO_DIR"

            rm -f \
                "$PYRO_GFF" \
                "$PYRO_NT" \
                "$PYRO_AA" \
                "$PYRO_QC" \
                "$PYRO_FLAG"

            cat > "$PYRO_JOB" <<EOL
#!/bin/bash

#SBATCH --cpus-per-task=${PYRO_CPUS}
#SBATCH --mem=${PYRO_MEM}
#SBATCH --job-name=Pyrodigal_${sample_name}
#SBATCH --error=${JOB_DIR}/${sample_name}_Pyrodigal_%j.error
#SBATCH --output=${JOB_DIR}/${sample_name}_Pyrodigal_%j.log
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=${PYRO_TIME}
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module load apps/anaconda3/2024.06/bin

source activate ${PYRODIGAL_ENV}

CONTIGS="${contigs_file}"

OUT_DIR="${PYRO_DIR}"

GFF="${PYRO_GFF}"

NT="${PYRO_NT}"

AA="${PYRO_AA}"

ASSEMBLY_QC="${ASSEMBLY_QC}"

PREV_FLAG="${metaspades_flag}"

SUCCESS_FLAG="${PYRO_FLAG}"

PYRO_QC="${PYRO_QC}"

# =========================================================================
# Failure handling
# =========================================================================

cleanup_on_failure() {

    rc=$?

    if [[ $rc -ne 0 ]]; then
        echo
        echo "ERROR: eggNOG-mapper pipeline failed for ${sample_name}."
        echo "The failure may have occurred during eggNOG execution or output validation."
        rm -f "$SUCCESS_FLAG"
    fi

    exit $rc
}

trap cleanup_on_failure EXIT

# =========================================================================
# Verify MetaSPAdes
# =========================================================================

if [[ ! -f "\$PREV_FLAG" ]]; then

    echo "ERROR: MetaSPAdes completion flag missing:"
    echo "\$PREV_FLAG"

    exit 1

fi

# =========================================================================
# Verify assembly
# =========================================================================

if [[ ! -s "\$CONTIGS" ]]; then

    echo "ERROR: Assembly missing or empty:"
    echo "\$CONTIGS"

    exit 1

fi

if [[ ! -s "\$ASSEMBLY_QC" ]]; then

    echo "ERROR: Assembly QC file missing:"
    echo "\$ASSEMBLY_QC"

    exit 1

fi

# =========================================================================
# Run Pyrodigal
# =========================================================================

echo
echo "============================================================"
echo "Pyrodigal"
echo "============================================================"

echo "Sample:"
echo "  ${sample_name}"

echo "Assembly:"
echo "  \$CONTIGS"

echo "Output:"
echo "  \$OUT_DIR"

echo "Threads:"
echo "  \$SLURM_CPUS_PER_TASK"

echo

mkdir -p "\$OUT_DIR"

pyrodigal \
    -i "\$CONTIGS" \
    -o "\$GFF" \
    -a "\$AA" \
    -d "\$NT" \
    -p meta \
    -j "\$SLURM_CPUS_PER_TASK"

# =========================================================================
# Validate Pyrodigal outputs
# =========================================================================

for file in "\$GFF" "\$NT" "\$AA"; do

    if [[ ! -s "\$file" ]]; then

        echo "ERROR: Missing/empty Pyrodigal output:"
        echo "  \$file"

        exit 1

    fi

done

# =========================================================================
# Pyrodigal output consistency QC
# =========================================================================

python - "\$GFF" "\$NT" "\$AA" "\$PYRO_QC" "${sample_name}" <<'PY'

import sys

gff = sys.argv[1]
nt = sys.argv[2]
aa = sys.argv[3]
qc = sys.argv[4]
sample = sys.argv[5]

def fasta_stats(path):

    count = 0
    total = 0
    ids = set()
    duplicate = 0
    current = None
    seq = []

    with open(path) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if current is not None:

                    sequence = "".join(seq)

                    if not sequence:
                        raise RuntimeError(
                            f"Empty sequence in {path}: {current}"
                        )

                    count += 1
                    total += len(sequence)

                current = line[1:].split()[0]
                seq = []

                if current in ids:
                    duplicate += 1

                ids.add(current)

            else:

                if current is None:
                    raise RuntimeError(
                        f"Sequence before FASTA header in {path}"
                    )

                seq.append(line)

    if current is not None:

        sequence = "".join(seq)

        if not sequence:
            raise RuntimeError(
                f"Empty sequence in {path}: {current}"
            )

        count += 1
        total += len(sequence)

    return count, total, duplicate


nt_count, nt_bp, nt_dup = fasta_stats(nt)

aa_count, aa_aa, aa_dup = fasta_stats(aa)

# -------------------------------------------------------------------------
# Count CDS features in GFF
# -------------------------------------------------------------------------

gff_cds = 0

with open(gff) as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) >= 3 and fields[2] == "CDS":
            gff_cds += 1


# -------------------------------------------------------------------------
# Validate consistency
# -------------------------------------------------------------------------

errors = []

if nt_count == 0:
    errors.append("No nucleotide predictions found.")

if aa_count == 0:
    errors.append("No protein predictions found.")

if gff_cds == 0:
    errors.append("No CDS features found in GFF.")

if nt_count != aa_count:
    errors.append(
        f"NT/protein prediction count mismatch: "
        f"{nt_count} vs {aa_count}"
    )

if gff_cds != aa_count:
    errors.append(
        f"GFF CDS/protein count mismatch: "
        f"{gff_cds} vs {aa_count}"
    )

if nt_dup > 0:
    errors.append(
        f"{nt_dup} duplicate nucleotide IDs."
    )

if aa_dup > 0:
    errors.append(
        f"{aa_dup} duplicate protein IDs."
    )

if errors:

    print(
        "ERROR: Pyrodigal validation failed:",
        file=sys.stderr
    )

    for error in errors:

        print(
            f"  - {error}",
            file=sys.stderr
        )

    sys.exit(1)

# -------------------------------------------------------------------------
# Write QC
# -------------------------------------------------------------------------

with open(qc, "w") as out:

    out.write(
        "Sample\tCDS\tNucleotide_sequences\tProtein_sequences\t"
        "Predicted_nt_bp\tPredicted_aa_residues\n"
    )

    out.write(
        f"{sample}\t"
        f"{gff_cds}\t"
        f"{nt_count}\t"
        f"{aa_count}\t"
        f"{nt_bp}\t"
        f"{aa_aa}\n"
    )

print("Pyrodigal QC passed.")
print(f"  CDS:       {gff_cds:,}")
print(f"  Proteins:  {aa_count:,}")
print(f"  NT seqs:   {nt_count:,}")

PY

# =========================================================================
# Final validation
# =========================================================================

if [[ ! -s "\$PYRO_QC" ]]; then

    echo "ERROR: Pyrodigal QC file missing."

    exit 1

fi

# =========================================================================
# Success flag
# =========================================================================

printf "%s\n" "\$(date --iso-8601=seconds)" > "\$SUCCESS_FLAG"

echo
echo "============================================================"
echo "Pyrodigal completed successfully"
echo "============================================================"

echo "Success flag:"
echo "  \$SUCCESS_FLAG"

EOL

            chmod +x "$PYRO_JOB"

            # -------------------------------------------------------------
            # Submit Pyrodigal
            # -------------------------------------------------------------

            echo
            echo "Submitting Pyrodigal:"
            echo "  $PYRO_JOB"

            PYRO_SUBMIT_OUTPUT=$(sbatch "$PYRO_JOB")

            echo "$PYRO_SUBMIT_OUTPUT"

            PYRO_JOB_ID=$(echo "$PYRO_SUBMIT_OUTPUT" | awk '{print $NF}')

            if [[ ! "$PYRO_JOB_ID" =~ ^[0-9]+$ ]]; then

                echo
                echo "ERROR: Could not determine Pyrodigal SLURM job ID."

                exit 1

            fi

            echo
            echo "Pyrodigal job ID:"
            echo "  $PYRO_JOB_ID"

        else

            echo
            echo "Pyrodigal already complete."

            PYRO_JOB_ID=""

        fi

        # =================================================================
        # 5J. Generate eggNOG job
        # =================================================================

        if [[ "$EGGNOG_COMPLETE" -eq 0 ]]; then

            echo
            echo "Preparing eggNOG-mapper job."

            mkdir -p "$EGGNOG_DIR"

            # -------------------------------------------------------------
            # Do NOT remove existing eggNOG results here if a Pyrodigal
            # job is still running. The eggNOG job will validate its input.
            # -------------------------------------------------------------

            rm -f "$EGGNOG_FLAG"

            # -------------------------------------------------------------
            # Determine dependency.
            #
            # If Pyrodigal was just submitted, wait for it.
            #
            # If Pyrodigal was already complete, eggNOG can be submitted
            # immediately.
            # -------------------------------------------------------------

            DEPENDENCY=""

            if [[ -n "$PYRO_JOB_ID" ]]; then

                DEPENDENCY="--dependency=afterok:${PYRO_JOB_ID}"

            fi

            cat > "$EGGNOG_JOB" <<EOL
#!/bin/bash

#SBATCH --cpus-per-task=${EGGNOG_CPUS}
#SBATCH --mem=${EGGNOG_MEM}
#SBATCH --job-name=EggnogMapper_${sample_name}
#SBATCH --error=${JOB_DIR}/${sample_name}_EggNOG_%j.error
#SBATCH --output=${JOB_DIR}/${sample_name}_EggNOG_%j.log
#SBATCH --partition=k2-medpri,k2-bioinf,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=${EGGNOG_TIME}
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module load apps/anaconda3/2024.06/bin

source activate ${EGGNOG_ENV}

export EGGNOG_DATA_DIR="${EGGNOG_DB}"

AA="${PYRO_AA}"

PYRO_FLAG="${PYRO_FLAG}"

PYRO_QC="${PYRO_QC}"

OUT_DIR="${EGGNOG_DIR}"

PREFIX="${EGGNOG_PREFIX}"

XLSX="${EGGNOG_XLSX}"

ANNOTATIONS="${EGGNOG_ANNOTATIONS}"

SUCCESS_FLAG="${EGGNOG_FLAG}"

QC="${EGGNOG_QC}"

# =========================================================================
# Failure handling
# =========================================================================

cleanup_on_failure() {

    rc=\$?

    if [[ \$rc -ne 0 ]]; then

        echo
        echo "ERROR: eggNOG-mapper failed for ${sample_name}."

        rm -f "\$SUCCESS_FLAG"

    fi

    exit \$rc
}

trap cleanup_on_failure EXIT

# =========================================================================
# Verify Pyrodigal completion
# =========================================================================

if [[ ! -f "\$PYRO_FLAG" ]]; then

    echo "ERROR: Pyrodigal success flag not found:"
    echo "\$PYRO_FLAG"

    exit 1

fi

# =========================================================================
# Verify Pyrodigal QC
# =========================================================================

if [[ ! -s "\$PYRO_QC" ]]; then

    echo "ERROR: Pyrodigal QC file missing:"
    echo "\$PYRO_QC"

    exit 1

fi

# =========================================================================
# Verify protein input
# =========================================================================

if [[ ! -s "\$AA" ]]; then

    echo "ERROR: Pyrodigal protein FASTA missing or empty:"
    echo "\$AA"

    exit 1

fi

# =========================================================================
# Verify that protein FASTA contains sequences
# =========================================================================

PROTEIN_COUNT=\$(grep -c '^>' "\$AA" || true)

if [[ "\$PROTEIN_COUNT" -eq 0 ]]; then

    echo "ERROR: No protein sequences found in:"
    echo "\$AA"

    exit 1

fi

echo
echo "============================================================"
echo "eggNOG-mapper"
echo "============================================================"

echo "Sample:"
echo "  ${sample_name}"

echo "Input protein FASTA:"
echo "  \$AA"

echo "Protein sequences:"
echo "  \$PROTEIN_COUNT"

echo "Output directory:"
echo "  \$OUT_DIR"

echo "CPUs:"
echo "  \$SLURM_CPUS_PER_TASK"

echo

mkdir -p "\$OUT_DIR"

# =========================================================================
# Run eggNOG-mapper
# =========================================================================

emapper.py \
    -i "\$AA" \
    --output_dir "\$OUT_DIR" \
    --output "\$PREFIX" \
    --score 60 \
    --subject_cover 60 \
    --sensmode sensitive \
    --dbmem \
    --decorate_gff yes \
    --excel \
    --cpu "\$SLURM_CPUS_PER_TASK" \
    --override

echo
echo "Checking eggNOG output files:"
echo

for file in "$ANNOTATIONS" "$XLSX"; do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Missing or empty eggNOG output:"
        echo "  $file"
        exit 1
    fi

    echo "  OK: $(ls -lh "$file")"
done

# =========================================================================
# Validate eggNOG outputs
# =========================================================================

if [[ ! -s "\$ANNOTATIONS" ]]; then

    echo "ERROR: eggNOG annotation output missing or empty:"
    echo "\$ANNOTATIONS"

    exit 1

fi

if [[ ! -s "\$XLSX" ]]; then

    echo "ERROR: eggNOG Excel output missing or empty:"
    echo "\$XLSX"

    exit 1

fi

# =========================================================================
# Basic annotation statistics
# =========================================================================

python - "\$ANNOTATIONS" "\$QC" "${sample_name}" "\$PROTEIN_COUNT" <<'PY'

import sys

annotations = sys.argv[1]
qc = sys.argv[2]
sample = sys.argv[3]
protein_count = int(sys.argv[4])

rows = 0
annotated = 0
unannotated = 0

with open(annotations) as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 1:
            continue

        rows += 1

        # eggNOG-mapper uses "-" for unavailable annotations
        #
        # We count a row as annotated if it contains a non-placeholder
        # annotation in any of the major functional fields.

        useful = False

        for field in fields[1:]:

            value = field.strip()

            if value and value != "-":

                useful = True
                break

        if useful:
            annotated += 1
        else:
            unannotated += 1


# -------------------------------------------------------------------------
# Write QC
# -------------------------------------------------------------------------

with open(qc, "w") as out:

    out.write(
        "Sample\tInput_proteins\tAnnotation_rows\t"
        "Rows_with_annotation\tRows_without_annotation\t"
        "Percent_with_annotation\n"
    )

    percentage = (
        annotated / rows * 100
        if rows > 0
        else 0
    )

    out.write(
        f"{sample}\t"
        f"{protein_count}\t"
        f"{rows}\t"
        f"{annotated}\t"
        f"{unannotated}\t"
        f"{percentage:.4f}\n"
    )

print("eggNOG annotation QC:")
print(f"  Input proteins:       {protein_count:,}")
print(f"  Annotation rows:      {rows:,}")
print(f"  Rows with annotation: {annotated:,}")
print(f"  Rows without:         {unannotated:,}")
print(f"  Annotated:            {percentage:.2f}%")

# -------------------------------------------------------------------------
# Do not fail solely because some proteins lack annotations.
# -------------------------------------------------------------------------

if rows == 0:

    print(
        "ERROR: eggNOG produced zero annotation rows.",
        file=sys.stderr
    )

    sys.exit(1)

PY

# =========================================================================
# Final eggNOG success flag
# =========================================================================

printf "%s\n" "\$(date --iso-8601=seconds)" > "\$SUCCESS_FLAG"

echo
echo "============================================================"
echo "eggNOG-mapper completed successfully"
echo "============================================================"

echo "Success flag:"
echo "  \$SUCCESS_FLAG"

EOL

            chmod +x "$EGGNOG_JOB"

            # -------------------------------------------------------------
            # Submit eggNOG
            # -------------------------------------------------------------

            echo
            echo "Submitting eggNOG-mapper:"
            echo "  $EGGNOG_JOB"

            if [[ -n "$PYRO_JOB_ID" ]]; then

                echo
                echo "Dependency:"
                echo "  afterok:$PYRO_JOB_ID"

                EGGNOG_SUBMIT_OUTPUT=$(sbatch \
                    --dependency="afterok:${PYRO_JOB_ID}" \
                    "$EGGNOG_JOB")

            else

                EGGNOG_SUBMIT_OUTPUT=$(sbatch \
                    "$EGGNOG_JOB")

            fi

            echo "$EGGNOG_SUBMIT_OUTPUT"

        else

            echo
            echo "eggNOG-mapper already complete."

        fi

    done

done

# =========================================================================
# 6. Final message
# =========================================================================

echo
echo "============================================================"
echo "Preparation/submission complete"
echo "============================================================"

echo
echo "The pipeline uses:"
echo
echo "  MetaSPAdes"
echo "       |"
echo "       v"
echo "  Assembly QC"
echo "       |"
echo "       v"
echo "  Pyrodigal"
echo "       |"
echo "       v"
echo "  Pyrodigal QC"
echo "       |"
echo "       v"
echo "  eggNOG-mapper"
echo "       |"
echo "       v"
echo "  eggNOG QC"
echo
echo "Existing successful samples were skipped."
echo
echo "============================================================"
echo "Done."
echo "============================================================"
c