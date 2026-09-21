#!/bin/bash
#SBATCH --job-name=Prepare_ReadMap_MethanogenGenomes_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_ReadMap_MethanogenGenomes_Jobs_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_ReadMap_MethanogenGenomes_Jobs_Ran_%j.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
shopt -s nullglob

module load apps/anaconda3/2024.06/bin

WORK_DIR="/mnt/scratch2/users/3057556/Ran/40-1310848928/00_fastq"
JOB_DIR="/users/3057556/jobs/readmap_target_genomes_jobs_ran"
OUTDIR_NAME="27_Methanogens_Readmapped"

# Pre-built Bowtie2 index containing the 27 methanogen reference genomes.
BT_DB="/mnt/scratch2/igfs-databases/27_methanogen_genomes/bt2_db/all_27_methanogens"

mkdir -p "$JOB_DIR"
cd "$WORK_DIR" || exit 1

# Validate the pre-built Bowtie2 index before generating jobs.
if [[ ! -s "${BT_DB}.1.bt2" && ! -s "${BT_DB}.1.bt2l" ]]; then
    echo "ERROR: Bowtie2 index not found for $BT_DB" >&2
    echo "Expected ${BT_DB}.1.bt2 or ${BT_DB}.1.bt2l" >&2
    exit 1
fi

for dir in ./*; do
    [[ -d "$dir" ]] || continue

    r1_files=("$dir"/*_trimmed_unpaired.fastq.gz*)
    [[ ${#r1_files[@]} -gt 0 ]] || continue

    for reads_1 in "${r1_files[@]}"; do
        filename=$(basename "$reads_1")
        sample_name=$(basename "$dir")

        # Replace only the paired-end marker.
        reads_2_name="${filename/_R1/_R2}"
        reads_2="$dir/$reads_2_name"

        if [[ ! -s "$reads_2" ]]; then
            echo "WARNING: R2 missing for $sample_name: $reads_2"
            continue
        fi

        outdir="$WORK_DIR/$OUTDIR_NAME/$sample_name"
        mkdir -p "$outdir"

        bam="$outdir/${sample_name}_readmapped_sorted.bam"
        bai="${bam}.bai"
        flagstat="$outdir/${sample_name}_readmapped_flagstat.txt"
        idxstats="$outdir/${sample_name}_readmapped_idxstats.txt"
        coverage="$outdir/${sample_name}_reference_coverage.tsv"
        mapping_summary="$outdir/${sample_name}_mapping_summary.tsv"
        bowtie_log="$outdir/${sample_name}_bowtie2.log"
        done_flag="$outdir/${sample_name}_readmapped_done.flag"

        job_script="$JOB_DIR/${sample_name}_readmap_job.sh"

        if [[ -s "$done_flag" &&
              -s "$bam" &&
              -s "$bai" &&
              -s "$flagstat" &&
              -s "$idxstats" &&
              -s "$coverage" &&
              -s "$mapping_summary" ]]; then
            echo "Skipping $sample_name - completed mapping already exists."
            continue
        fi

        echo "============================================================"
        echo "Sample:            $sample_name"
        echo "R1:                $reads_1"
        echo "R2:                $reads_2"
        echo "Output directory:  $outdir"
        echo "Bowtie2 database:   $BT_DB"
        echo "Sorted BAM:        $bam"
        echo "Job script:        $job_script"
        echo "============================================================"

        cat <<EOL > "$job_script"
#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=ReadMapMethanogen_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

THREADS="\${SLURM_CPUS_PER_TASK}"

SAMPLE="$sample_name"
READS_1="$dir/$filename"
READS_2="$dir/$reads_2_name"
BT_DB="$BT_DB"

OUTDIR="$outdir"
BAM="$bam"
BAI="$bai"
FLAGSTAT="$flagstat"
IDXSTATS="$idxstats"
COVERAGE="$coverage"
MAPPING_SUMMARY="$mapping_summary"
BOWTIE_LOG="$bowtie_log"
DONE_FLAG="$done_flag"

TMP_BAM="\$OUTDIR/\${SAMPLE}_readmapped_sorted.bam.tmp"

mkdir -p "\$OUTDIR"
rm -f "\$DONE_FLAG"

for cmd in bowtie2 samtools awk; do
    command -v "\$cmd" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found: \$cmd" >&2
        exit 1
    }
done

[[ -s "\$READS_1" ]] || { echo "ERROR: R1 missing: \$READS_1" >&2; exit 1; }
[[ -s "\$READS_2" ]] || { echo "ERROR: R2 missing: \$READS_2" >&2; exit 1; }
[[ -s "\${BT_DB}.1.bt2" || -s "\${BT_DB}.1.bt2l" ]] || {
    echo "ERROR: Bowtie2 index missing: \$BT_DB" >&2
    exit 1
}

# Map reads. The temporary BAM is only renamed into place after sorting
# completes successfully, preventing an incomplete BAM looking valid.
if [[ ! -s "\$BAM" || ! -s "\$BAI" ]]; then
    rm -f "\$TMP_BAM" "\$BAM" "\$BAI"

    echo "Starting Bowtie2 mapping for \$SAMPLE..."
    date

    bowtie2 \
        -x "\$BT_DB" \
        -1 "\$READS_1" \
        -2 "\$READS_2" \
        --very-sensitive-local \
        -p "\$THREADS" \
        2> "\$BOWTIE_LOG" \
    | samtools view -@ "\$THREADS" -b - \
    | samtools sort -@ "\$THREADS" -o "\$TMP_BAM" -

    mv "\$TMP_BAM" "\$BAM"
    samtools index -@ "\$THREADS" "\$BAM"
fi

[[ -s "\$BAM" ]] || { echo "ERROR: BAM missing." >&2; exit 1; }
[[ -s "\$BAI" ]] || { echo "ERROR: BAM index missing." >&2; exit 1; }

# Mapping QC.
samtools flagstat "\$BAM" > "\$FLAGSTAT"
samtools idxstats "\$BAM" > "\$IDXSTATS"

# Reference-level coverage. If a methanogen genome consists of multiple
# reference records, these remain separate here and can later be aggregated
# to genome/accession level by MetaPont.
samtools coverage "\$BAM" \
    | awk '
        BEGIN {
            OFS="\t"
            print "Reference","Length","Mapped_Reads","Covered_Bases",
                  "Breadth","Mean_Depth","Mean_BaseQ","Mean_MapQ"
        }
        NR > 1 {
            print \$1,\$3-\$2+1,\$4,\$5,\$6,\$7,\$8,\$9
        }
    ' > "\$COVERAGE"

# Compact sample-level mapping summary.
TOTAL_PRIMARY_ALIGNMENTS=\$(samtools view -c -F 2304 "\$BAM")
MAPPED_PRIMARY_ALIGNMENTS=\$(samtools view -c -F 2308 "\$BAM")

if [[ "\$TOTAL_PRIMARY_ALIGNMENTS" -gt 0 ]]; then
    MAPPED_PROPORTION=\$(awk -v m="\$MAPPED_PRIMARY_ALIGNMENTS" -v t="\$TOTAL_PRIMARY_ALIGNMENTS" \
        'BEGIN {printf "%.8f", m/t}')
else
    MAPPED_PROPORTION="0"
fi

REFERENCES_TOTAL=\$(awk '\$1 != "*" {n++} END {print n+0}' "\$IDXSTATS")
REFERENCES_SUPPORTED=\$(awk 'NR > 1 && \$3 > 0 {n++} END {print n+0}' "\$COVERAGE")
IDXSTATS_MAPPED=\$(awk '\$1 != "*" {sum += \$3} END {print sum+0}' "\$IDXSTATS")

cat > "\$MAPPING_SUMMARY" <<SUMMARY
Field	Value
Sample	\$SAMPLE
Bowtie2_database	\$BT_DB
Mapping_mode	--very-sensitive-local
Threads	\$THREADS
Total_primary_alignments	\$TOTAL_PRIMARY_ALIGNMENTS
Mapped_primary_alignments	\$MAPPED_PRIMARY_ALIGNMENTS
Mapped_alignment_proportion	\$MAPPED_PROPORTION
Mapped_alignments_idxstats	\$IDXSTATS_MAPPED
Reference_sequences_total	\$REFERENCES_TOTAL
Reference_sequences_supported	\$REFERENCES_SUPPORTED
SUMMARY

# Final validation.
[[ -s "\$FLAGSTAT" ]] || { echo "ERROR: flagstat is empty." >&2; exit 1; }
[[ -s "\$IDXSTATS" ]] || { echo "ERROR: idxstats is empty." >&2; exit 1; }
[[ -s "\$COVERAGE" ]] || { echo "ERROR: coverage table is empty." >&2; exit 1; }
[[ -s "\$MAPPING_SUMMARY" ]] || { echo "ERROR: mapping summary is empty." >&2; exit 1; }

printf "Completed %s\n" "\$(date -Is)" > "\$DONE_FLAG"

echo "============================================================"
echo "Mapping complete:      \$SAMPLE"
echo "Primary alignments:    \$TOTAL_PRIMARY_ALIGNMENTS"
echo "Mapped alignments:     \$MAPPED_PRIMARY_ALIGNMENTS"
echo "Mapping proportion:    \$MAPPED_PROPORTION"
echo "References supported:  \$REFERENCES_SUPPORTED / \$REFERENCES_TOTAL"
echo "============================================================"
EOL

        chmod 750 "$job_script"
        echo "Submitting job for $sample_name..."
        sbatch "$job_script"
    done
done
