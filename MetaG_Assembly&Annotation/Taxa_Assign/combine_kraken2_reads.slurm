#!/bin/bash

#SBATCH --job-name=Combine_Kraken2_Reads_Ran
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --error=/users/3057556/jobs/Combine_Kraken2_Reads_Ran_%j.error
#SBATCH --output=/users/3057556/jobs/Combine_Kraken2_Reads_Ran_%j.log
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# =========================================================================
# Kraken2 aggregation, taxonomic summaries and visualisations
#
# Designed for ~80 samples.
#
# IMPORTANT:
# Kraken2 uses abbreviated rank codes:
#
#   K = kingdom
#   P = phylum
#   C = class
#   O = order
#   F = family
#   G = genus
#   S = species
#
# Ranks such as G1, G2, F1 etc. are subordinate Kraken2 ranks and are
# deliberately excluded from the primary rank-level abundance tables.
#
# The script is restartable:
#   - Existing valid outputs are not regenerated.
#   - Existing combined Kraken2 results are reused.
#   - Individual visualisations are generated independently.
#   - Krona is generated per sample rather than as one enormous HTML.
# =========================================================================

set -euo pipefail

# =========================================================================
# 1. Load environment
# =========================================================================

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/kraken2_2.1.3

# =========================================================================
# 2. Define paths
# =========================================================================

WORK_DIR="/mnt/scratch2/users/3057556/Ran/40-1310848928/00_fastq"

OUTPUT_DIR="${WORK_DIR}/Kraken2_reads_combined"

COMBINED_DIR="${OUTPUT_DIR}/combined"
SUMMARY_DIR="${OUTPUT_DIR}/summaries"
ABUNDANCE_DIR="${OUTPUT_DIR}/abundance"
VIS_DIR="${OUTPUT_DIR}/visualisations"
KRONA_DIR="${OUTPUT_DIR}/Krona"

mkdir -p \
    "$OUTPUT_DIR" \
    "$COMBINED_DIR" \
    "$SUMMARY_DIR" \
    "$ABUNDANCE_DIR" \
    "$VIS_DIR" \
    "$KRONA_DIR"

cd "$WORK_DIR"

echo
echo "============================================================"
echo "Kraken2 aggregation"
echo "============================================================"
echo
echo "Working directory:"
echo "  $WORK_DIR"
echo
echo "Output directory:"
echo "  $OUTPUT_DIR"
echo

# =========================================================================
# 3. Locate Kraken2 files
# =========================================================================

mapfile -t KRAKEN_FILES < <(
    find "$WORK_DIR" \
        -type f \
        -name "*_kraken2_reads.txt" \
        ! -path "$OUTPUT_DIR/*" \
        | sort
)

if [[ "${#KRAKEN_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No Kraken2 read-level files found."
    exit 1
fi

mapfile -t REPORT_FILES < <(
    find "$WORK_DIR" \
        -type f \
        -name "*_kraken2_reads_report.txt" \
        ! -path "$OUTPUT_DIR/*" \
        | sort
)

if [[ "${#REPORT_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No Kraken2 report files found."
    exit 1
fi

echo "Found ${#KRAKEN_FILES[@]} read-level files."
echo "Found ${#REPORT_FILES[@]} report files."

# =========================================================================
# 4. Define combined outputs
# =========================================================================

COMBINED_READS="${COMBINED_DIR}/Ran_Kraken2_combined_read_level.tsv.gz"

COMBINED_REPORT="${COMBINED_DIR}/Ran_Kraken2_combined_report.tsv.gz"

SAMPLE_SUMMARY="${SUMMARY_DIR}/Ran_Kraken2_sample_summary.tsv"

# =========================================================================
# 5. Combine read-level Kraken2 results
# =========================================================================

if [[ -s "$COMBINED_READS" ]]; then

    echo
    echo "[SKIP] Combined read-level results already exist."

else

    echo
    echo "============================================================"
    echo "Combining read-level Kraken2 results"
    echo "============================================================"

    TMP_READS="${COMBINED_DIR}/Ran_Kraken2_combined_read_level.tsv.tmp"

    rm -f "$TMP_READS"

    printf "Sample\tKraken_classification\tTaxid\tRead_ID\tLength\tLCA_mapping\n" \
        > "$TMP_READS"

    for file in "${KRAKEN_FILES[@]}"; do

        sample_name=$(basename "$file" "_kraken2_reads.txt")

        echo "  Adding: $sample_name"

        awk -v sample="$sample_name" '
        BEGIN {
            OFS="\t"
        }

        {
            print sample, $1, $3, $2, $4, $5
        }
        ' "$file" >> "$TMP_READS"

    done

    gzip -c "$TMP_READS" > "$COMBINED_READS"

    rm -f "$TMP_READS"

    echo "Created:"
    echo "  $COMBINED_READS"

fi

# =========================================================================
# 6. Combine Kraken2 reports
# =========================================================================

if [[ -s "$COMBINED_REPORT" ]]; then

    echo
    echo "[SKIP] Combined Kraken2 report already exists."

else

    echo
    echo "============================================================"
    echo "Combining Kraken2 reports"
    echo "============================================================"

    TMP_REPORT="${COMBINED_DIR}/Ran_Kraken2_combined_report.tsv.tmp"

    rm -f "$TMP_REPORT"

    printf "Sample\tPercentage\tClade_reads\tTaxon_reads\tRank\tTaxID\tName\n" \
        > "$TMP_REPORT"

    for file in "${REPORT_FILES[@]}"; do

        sample_name=$(basename "$file" "_kraken2_reads_report.txt")

        echo "  Adding: $sample_name"

        awk -v sample="$sample_name" '
        BEGIN {
            FS="\t"
            OFS="\t"
        }

        NF >= 6 {

            percentage=$1
            clade_reads=$2
            taxon_reads=$3
            rank=$4
            taxid=$5
            name=$6

            sub(/^[ \t]+/, "", name)

            print sample,
                  percentage,
                  clade_reads,
                  taxon_reads,
                  rank,
                  taxid,
                  name
        }
        ' "$file" >> "$TMP_REPORT"

    done

    gzip -c "$TMP_REPORT" > "$COMBINED_REPORT"

    rm -f "$TMP_REPORT"

    echo "Created:"
    echo "  $COMBINED_REPORT"

fi

# =========================================================================
# 7. Sample-level summary
# =========================================================================

if [[ -s "$SAMPLE_SUMMARY" ]]; then

    echo
    echo "[SKIP] Sample summary already exists."

else

    echo
    echo "============================================================"
    echo "Generating sample summary"
    echo "============================================================"

    TMP_SUMMARY="${SUMMARY_DIR}/Ran_Kraken2_sample_summary.tsv.tmp"

    rm -f "$TMP_SUMMARY"

    printf "Sample\tTotal_reads\tClassified_reads\tUnclassified_reads\tPercent_classified\tPercent_unclassified\n" \
        > "$TMP_SUMMARY"

    for file in "${KRAKEN_FILES[@]}"; do

        sample_name=$(basename "$file" "_kraken2_reads.txt")

        awk -v sample="$sample_name" '
        BEGIN {
            total=0
            classified=0
            unclassified=0
        }

        {
            total++

            if ($1 == "C")
                classified++

            else if ($1 == "U")
                unclassified++
        }

        END {

            if (total > 0) {

                pc=(classified/total)*100
                pu=(unclassified/total)*100

                printf "%s\t%d\t%d\t%d\t%.4f\t%.4f\n",
                    sample,
                    total,
                    classified,
                    unclassified,
                    pc,
                    pu
            }
        }
        ' "$file" >> "$TMP_SUMMARY"

    done

    {
        head -n 1 "$TMP_SUMMARY"
        tail -n +2 "$TMP_SUMMARY" | sort -k1,1
    } > "${TMP_SUMMARY}.sorted"

    mv "${TMP_SUMMARY}.sorted" "$SAMPLE_SUMMARY"

    rm -f "$TMP_SUMMARY"

    echo "Created:"
    echo "  $SAMPLE_SUMMARY"

fi

# =========================================================================
# 8. Generate rank-specific abundance tables
# =========================================================================

echo
echo "============================================================"
echo "Generating taxonomic abundance tables"
echo "============================================================"

python - "$COMBINED_REPORT" "$SAMPLE_SUMMARY" "$ABUNDANCE_DIR" <<'PY'

import sys
import os
import gzip
from collections import defaultdict

report_file = sys.argv[1]
summary_file = sys.argv[2]
output_root = sys.argv[3]

# -------------------------------------------------------------------------
# Kraken2 primary taxonomic ranks
#
# Only the exact rank codes below are retained.
# G1/G2/etc. are NOT treated as genera.
# -------------------------------------------------------------------------

RANKS = {
    "K": "kingdom",
    "P": "phylum",
    "C": "class",
    "O": "order",
    "F": "family",
    "G": "genus",
    "S": "species",
}

# -------------------------------------------------------------------------
# Read total reads per sample
# -------------------------------------------------------------------------

sample_totals = {}

with open(summary_file) as f:

    next(f)

    for line in f:

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 2:
            continue

        sample = fields[0]

        try:
            total = int(fields[1])
        except ValueError:
            continue

        sample_totals[sample] = total


samples = sorted(sample_totals)


# -------------------------------------------------------------------------
# Store:
#
# counts[rank][taxid][sample] = clade reads
# names[rank][taxid] = taxon name
# -------------------------------------------------------------------------

counts = {
    rank: defaultdict(lambda: defaultdict(int))
    for rank in RANKS.values()
}

names = {
    rank: {}
    for rank in RANKS.values()
}


# -------------------------------------------------------------------------
# Read combined Kraken2 report
# -------------------------------------------------------------------------

with gzip.open(report_file, "rt") as f:

    next(f)

    for line in f:

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 7:
            continue

        sample = fields[0]
        clade_reads = fields[2]
        rank_code = fields[4]
        taxid = fields[5]
        name = fields[6]

        # Only primary ranks
        if rank_code not in RANKS:
            continue

        rank = RANKS[rank_code]

        try:
            clade_reads = int(clade_reads)
        except ValueError:
            continue

        name = name.strip()

        if not name:
            continue

        counts[rank][taxid][sample] = clade_reads
        names[rank][taxid] = name


# -------------------------------------------------------------------------
# Generate tables
# -------------------------------------------------------------------------

for rank in RANKS.values():

    rank_dir = os.path.join(output_root, rank)

    os.makedirs(rank_dir, exist_ok=True)

    count_file = os.path.join(
        rank_dir,
        f"Ran_Kraken2_{rank}_counts.tsv"
    )

    abundance_file = os.path.join(
        rank_dir,
        f"Ran_Kraken2_{rank}_relative_abundance.tsv"
    )

    presence_file = os.path.join(
        rank_dir,
        f"Ran_Kraken2_{rank}_presence.tsv"
    )

    metadata_file = os.path.join(
        rank_dir,
        f"Ran_Kraken2_{rank}_taxa_metadata.tsv"
    )

    # -------------------------------------------------------------
    # Sort taxa by total reads across samples
    # -------------------------------------------------------------

    taxon_totals = {}

    for taxid, sample_values in counts[rank].items():

        taxon_totals[taxid] = sum(
            sample_values.values()
        )

    ordered_taxids = sorted(
        taxon_totals,
        key=lambda taxid: (
            -taxon_totals[taxid],
            names[rank].get(taxid, "")
        )
    )

    # -------------------------------------------------------------
    # Counts
    # -------------------------------------------------------------

    with open(count_file, "w") as out:

        out.write(
            "TaxID\tTaxon\t" +
            "\t".join(samples) +
            "\tTotal\n"
        )

        for taxid in ordered_taxids:

            values = [
                counts[rank][taxid].get(sample, 0)
                for sample in samples
            ]

            out.write(
                f"{taxid}\t"
                f"{names[rank][taxid]}\t" +
                "\t".join(map(str, values)) +
                f"\t{sum(values)}\n"
            )

    # -------------------------------------------------------------
    # Relative abundance
    # -------------------------------------------------------------

    with open(abundance_file, "w") as out:

        out.write(
            "TaxID\tTaxon\t" +
            "\t".join(samples) +
            "\tMean_relative_abundance\n"
        )

        for taxid in ordered_taxids:

            values = []

            for sample in samples:

                total = sample_totals[sample]
                count = counts[rank][taxid].get(sample, 0)

                if total > 0:
                    value = (count / total) * 100.0
                else:
                    value = 0.0

                values.append(value)

            mean_value = (
                sum(values) / len(values)
                if values else 0.0
            )

            out.write(
                f"{taxid}\t"
                f"{names[rank][taxid]}\t" +
                "\t".join(
                    f"{x:.6f}" for x in values
                ) +
                f"\t{mean_value:.6f}\n"
            )

    # -------------------------------------------------------------
    # Presence / absence
    #
    # Any taxon with >0 clade reads is considered present.
    # -------------------------------------------------------------

    with open(presence_file, "w") as out:

        out.write(
            "TaxID\tTaxon\t" +
            "\t".join(samples) +
            "\n"
        )

        for taxid in ordered_taxids:

            values = [
                1 if counts[rank][taxid].get(sample, 0) > 0
                else 0
                for sample in samples
            ]

            out.write(
                f"{taxid}\t"
                f"{names[rank][taxid]}\t" +
                "\t".join(map(str, values)) +
                "\n"
            )

    # -------------------------------------------------------------
    # Metadata
    # -------------------------------------------------------------

    with open(metadata_file, "w") as out:

        out.write(
            "TaxID\tTaxon\tRank\tTotal_reads_across_samples\n"
        )

        for taxid in ordered_taxids:

            out.write(
                f"{taxid}\t"
                f"{names[rank][taxid]}\t"
                f"{rank}\t"
                f"{taxon_totals[taxid]}\n"
            )

    print(
        f"  {rank}: {len(ordered_taxids)} taxa"
    )

PY

# =========================================================================
# 9. Generate visualisations
# =========================================================================

echo
echo "============================================================"
echo "Generating visualisations"
echo "============================================================"

python - "$ABUNDANCE_DIR" "$SUMMARY_DIR" "$VIS_DIR" <<'PY'

import sys
import os
import csv

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


abundance_root = sys.argv[1]
summary_root = sys.argv[2]
vis_root = sys.argv[3]

os.makedirs(vis_root, exist_ok=True)


# =========================================================================
# Helpers
# =========================================================================

def read_tsv(path):

    with open(path, newline="") as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        rows = list(reader)

        return rows, reader.fieldnames


def save_if_missing(fig, path):

    if os.path.exists(path) and os.path.getsize(path) > 0:

        print(f"[SKIP] {path}")

    else:

        fig.savefig(
            path,
            dpi=200,
            bbox_inches="tight"
        )

        print(f"[CREATED] {path}")


# =========================================================================
# Sample summary
# =========================================================================

summary_file = os.path.join(
    summary_root,
    "Ran_Kraken2_sample_summary.tsv"
)

summary, _ = read_tsv(summary_file)

samples = [
    row["Sample"]
    for row in summary
]


# =========================================================================
# 1. Classification rate
# =========================================================================

outfile = os.path.join(
    vis_root,
    "classified_reads_per_sample.png"
)

if not os.path.exists(outfile):

    classified = [
        float(row["Percent_classified"])
        for row in summary
    ]

    unclassified = [
        float(row["Percent_unclassified"])
        for row in summary
    ]

    order = sorted(
        range(len(samples)),
        key=lambda i: classified[i]
    )

    plot_samples = [
        samples[i]
        for i in order
    ]

    plot_classified = [
        classified[i]
        for i in order
    ]

    plot_unclassified = [
        unclassified[i]
        for i in order
    ]

    fig, ax = plt.subplots(
        figsize=(16, 7)
    )

    ax.bar(
        range(len(plot_samples)),
        plot_classified,
        label="Classified"
    )

    ax.bar(
        range(len(plot_samples)),
        plot_unclassified,
        bottom=plot_classified,
        label="Unclassified"
    )

    ax.set_ylabel(
        "Percentage of reads"
    )

    ax.set_xlabel("Sample")

    ax.set_title(
        "Kraken2 classification rate"
    )

    ax.set_xticks(
        range(len(plot_samples))
    )

    ax.set_xticklabels(
        plot_samples,
        rotation=90,
        fontsize=7
    )

    ax.set_ylim(0, 100)

    ax.legend()

    fig.tight_layout()

    save_if_missing(
        fig,
        outfile
    )

    plt.close(fig)


# =========================================================================
# 2. Total reads
# =========================================================================

outfile = os.path.join(
    vis_root,
    "total_reads_per_sample.png"
)

if not os.path.exists(outfile):

    totals = [
        int(row["Total_reads"])
        for row in summary
    ]

    order = sorted(
        range(len(samples)),
        key=lambda i: totals[i]
    )

    plot_samples = [
        samples[i]
        for i in order
    ]

    plot_totals = [
        totals[i]
        for i in order
    ]

    fig, ax = plt.subplots(
        figsize=(16, 7)
    )

    ax.bar(
        range(len(plot_samples)),
        plot_totals
    )

    ax.set_ylabel(
        "Total reads"
    )

    ax.set_xlabel("Sample")

    ax.set_title(
        "Total reads per sample"
    )

    ax.set_xticks(
        range(len(plot_samples))
    )

    ax.set_xticklabels(
        plot_samples,
        rotation=90,
        fontsize=7
    )

    ax.set_yscale("log")

    fig.tight_layout()

    save_if_missing(
        fig,
        outfile
    )

    plt.close(fig)


# =========================================================================
# Function to generate stacked barplot
# =========================================================================

def make_stacked_barplot(
    rank,
    top_n,
    filename,
    title
):

    outfile = os.path.join(
        vis_root,
        filename
    )

    if os.path.exists(outfile):
        print(f"[SKIP] {outfile}")
        return

    abundance_file = os.path.join(
        abundance_root,
        rank,
        f"Ran_Kraken2_{rank}_relative_abundance.tsv"
    )

    if not os.path.exists(abundance_file):

        print(
            f"[WARNING] Missing {abundance_file}"
        )

        return

    rows, fields = read_tsv(
        abundance_file
    )

    if not rows:

        print(
            f"[WARNING] No taxa at {rank} level."
        )

        return

    sample_columns = [
        field
        for field in fields
        if field not in (
            "TaxID",
            "Taxon",
            "Mean_relative_abundance"
        )
    ]

    # Top taxa by mean abundance
    ranked = sorted(
        rows,
        key=lambda row: float(
            row["Mean_relative_abundance"]
        ),
        reverse=True
    )

    top_rows = ranked[:top_n]

    top_names = [
        row["Taxon"]
        for row in top_rows
    ]

    top_set = set(top_names)

    # Values for top taxa
    values = {}

    for row in top_rows:

        values[row["Taxon"]] = [
            float(row[sample])
            for sample in sample_columns
        ]

    # Everything else
    other = []

    for sample in sample_columns:

        value = 0.0

        for row in rows:

            if row["Taxon"] not in top_set:

                value += float(
                    row[sample]
                )

        other.append(value)

    fig, ax = plt.subplots(
        figsize=(18, 8)
    )

    bottom = [
        0.0
        for _ in sample_columns
    ]

    for taxon in top_names:

        vals = values[taxon]

        ax.bar(
            range(len(sample_columns)),
            vals,
            bottom=bottom,
            label=taxon
        )

        bottom = [
            x + y
            for x, y in zip(
                bottom,
                vals
            )
        ]

    ax.bar(
        range(len(sample_columns)),
        other,
        bottom=bottom,
        label="Other"
    )

    ax.set_ylabel(
        "Relative abundance (%)"
    )

    ax.set_xlabel("Sample")

    ax.set_title(title)

    ax.set_xticks(
        range(len(sample_columns))
    )

    ax.set_xticklabels(
        sample_columns,
        rotation=90,
        fontsize=7
    )

    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8
    )

    fig.tight_layout()

    save_if_missing(
        fig,
        outfile
    )

    plt.close(fig)


# =========================================================================
# 3. Top genera
# =========================================================================

make_stacked_barplot(
    rank="genus",
    top_n=20,
    filename="top_20_genera_stacked_barplot.png",
    title="Top 20 genera across samples"
)


# =========================================================================
# 4. Top families
# =========================================================================

make_stacked_barplot(
    rank="family",
    top_n=20,
    filename="top_20_families_stacked_barplot.png",
    title="Top 20 families across samples"
)


# =========================================================================
# 5. Top phyla
# =========================================================================

make_stacked_barplot(
    rank="phylum",
    top_n=15,
    filename="top_15_phyla_stacked_barplot.png",
    title="Top 15 phyla across samples"
)


# =========================================================================
# 6. Genus heatmap
# =========================================================================

outfile = os.path.join(
    vis_root,
    "top_30_genera_heatmap.png"
)

if not os.path.exists(outfile):

    genus_file = os.path.join(
        abundance_root,
        "genus",
        "Ran_Kraken2_genus_relative_abundance.tsv"
    )

    rows, fields = read_tsv(genus_file)

    if rows:

        sample_columns = [
            field
            for field in fields
            if field not in (
                "TaxID",
                "Taxon",
                "Mean_relative_abundance"
            )
        ]

        ranked = sorted(
            rows,
            key=lambda row: float(
                row["Mean_relative_abundance"]
            ),
            reverse=True
        )

        top_rows = ranked[:30]

        if top_rows:

            matrix = [
                [
                    float(row[sample])
                    for sample in sample_columns
                ]
                for row in top_rows
            ]

            names = [
                row["Taxon"]
                for row in top_rows
            ]

            fig_width = max(
                14,
                len(sample_columns) * 0.22
            )

            fig_height = max(
                8,
                len(names) * 0.30
            )

            fig, ax = plt.subplots(
                figsize=(
                    fig_width,
                    fig_height
                )
            )

            image = ax.imshow(
                matrix,
                aspect="auto",
                interpolation="nearest"
            )

            ax.set_yticks(
                range(len(names))
            )

            ax.set_yticklabels(
                names,
                fontsize=8
            )

            ax.set_xticks(
                range(len(sample_columns))
            )

            ax.set_xticklabels(
                sample_columns,
                rotation=90,
                fontsize=6
            )

            ax.set_xlabel("Sample")

            ax.set_ylabel("Genus")

            ax.set_title(
                "Top 30 genera by mean relative abundance"
            )

            fig.colorbar(
                image,
                ax=ax,
                label="Relative abundance (%)"
            )

            fig.tight_layout()

            save_if_missing(
                fig,
                outfile
            )

            plt.close(fig)

        else:

            print(
                "[WARNING] No genera available "
                "for heatmap."
            )

    else:

        print(
            "[WARNING] Genus table is empty."
        )


# =========================================================================
# 7. Genus richness
# =========================================================================

outfile = os.path.join(
    vis_root,
    "genus_richness_per_sample.png"
)

if not os.path.exists(outfile):

    genus_file = os.path.join(
        abundance_root,
        "genus",
        "Ran_Kraken2_genus_relative_abundance.tsv"
    )

    rows, fields = read_tsv(genus_file)

    sample_columns = [
        field
        for field in fields
        if field not in (
            "TaxID",
            "Taxon",
            "Mean_relative_abundance"
        )
    ]

    richness = []

    for sample in sample_columns:

        n = sum(
            1
            for row in rows
            if float(row[sample]) > 0
        )

        richness.append(n)

    order = sorted(
        range(len(sample_columns)),
        key=lambda i: richness[i]
    )

    plot_samples = [
        sample_columns[i]
        for i in order
    ]

    plot_richness = [
        richness[i]
        for i in order
    ]

    fig, ax = plt.subplots(
        figsize=(16, 7)
    )

    ax.bar(
        range(len(plot_samples)),
        plot_richness
    )

    ax.set_ylabel(
        "Number of detected genera"
    )

    ax.set_xlabel("Sample")

    ax.set_title(
        "Kraken2 genus-level richness"
    )

    ax.set_xticks(
        range(len(plot_samples))
    )

    ax.set_xticklabels(
        plot_samples,
        rotation=90,
        fontsize=7
    )

    fig.tight_layout()

    save_if_missing(
        fig,
        outfile
    )

    plt.close(fig)


print()
print("Visualisation generation complete.")

PY

# =========================================================================
# 10. Generate per-sample Krona
# =========================================================================

echo
echo "============================================================"
echo "Generating per-sample Krona files"
echo "============================================================"

if command -v ktImportText >/dev/null 2>&1; then

    echo "Krona detected:"
    echo "  $(command -v ktImportText)"
    echo

    for file in "${REPORT_FILES[@]}"; do

        sample_name=$(basename "$file" "_kraken2_reads_report.txt")

        KRONA_HTML="${KRONA_DIR}/${sample_name}_Krona.html"

        if [[ -s "$KRONA_HTML" ]]; then

            echo "[SKIP] Krona already exists: $sample_name"
            continue

        fi

        KRONA_INPUT="${KRONA_DIR}/${sample_name}_Krona_input.tsv"

        echo "Generating Krona input for: $sample_name"

        # -------------------------------------------------------------
        # Reconstruct Kraken2 hierarchy from indentation.
        #
        # First column = abundance/count.
        # Remaining columns = taxonomy hierarchy.
        #
        # We only need this for visualisation; the abundance tables
        # above are the authoritative summary outputs.
        # -------------------------------------------------------------

        awk '
        BEGIN {
            FS="\t"
            OFS="\t"
        }

        NF >= 6 {

            count=$2
            raw_name=$6

            # Number of leading spaces in taxon name
            match(raw_name, /[^ ]/)

            leading=RSTART-1

            depth=int(leading/2)

            name=raw_name
            sub(/^[ \t]+/, "", name)

            if (name == "")
                next

            names[depth]=name

            for (i=depth+1; i<100; i++)
                delete names[i]

            if (count <= 0)
                next

            taxonomy=""

            for (i=0; i<=depth; i++) {

                if (names[i] != "") {

                    if (taxonomy == "")
                        taxonomy=names[i]

                    else
                        taxonomy=taxonomy "\t" names[i]
                }
            }

            if (taxonomy != "")
                print count, taxonomy
        }

        ' "$file" > "$KRONA_INPUT"

        if [[ ! -s "$KRONA_INPUT" ]]; then

            echo "  WARNING: Empty Krona input."
            rm -f "$KRONA_INPUT"
            continue

        fi

        echo "  Generating HTML..."

        ktImportText \
            -o "$KRONA_HTML" \
            "$KRONA_INPUT"

        if [[ -s "$KRONA_HTML" ]]; then

            echo "  [OK] $KRONA_HTML"

        else

            echo "  WARNING: Krona HTML was not created."

        fi

    done

else

    echo
    echo "Krona not found."
    echo "Skipping Krona generation."
    echo
    echo "All non-Krona outputs will still be produced."

fi

# =========================================================================
# 11. Validation
# =========================================================================

echo
echo "============================================================"
echo "Validating outputs"
echo "============================================================"

FAILED=0


check_file() {

    local file="$1"

    if [[ -s "$file" ]]; then

        echo "[OK] $file"

    else

        echo "[FAIL] $file"

        FAILED=1

    fi
}


check_file "$COMBINED_READS"
check_file "$COMBINED_REPORT"
check_file "$SAMPLE_SUMMARY"


for rank in \
    kingdom \
    phylum \
    class \
    order \
    family \
    genus \
    species
do

    check_file \
        "${ABUNDANCE_DIR}/${rank}/Ran_Kraken2_${rank}_counts.tsv"

    check_file \
        "${ABUNDANCE_DIR}/${rank}/Ran_Kraken2_${rank}_relative_abundance.tsv"

    check_file \
        "${ABUNDANCE_DIR}/${rank}/Ran_Kraken2_${rank}_presence.tsv"

    check_file \
        "${ABUNDANCE_DIR}/${rank}/Ran_Kraken2_${rank}_taxa_metadata.tsv"

done


check_file \
    "${VIS_DIR}/classified_reads_per_sample.png"

check_file \
    "${VIS_DIR}/total_reads_per_sample.png"

check_file \
    "${VIS_DIR}/top_20_genera_stacked_barplot.png"

check_file \
    "${VIS_DIR}/top_20_families_stacked_barplot.png"

check_file \
    "${VIS_DIR}/top_15_phyla_stacked_barplot.png"

check_file \
    "${VIS_DIR}/top_30_genera_heatmap.png"

check_file \
    "${VIS_DIR}/genus_richness_per_sample.png"


# =========================================================================
# 12. Completion flag
# =========================================================================

SUCCESS_FLAG="${OUTPUT_DIR}/Ran_Kraken2_combination_done.flag"

if [[ "$FAILED" -eq 0 ]]; then

    printf "%s\n" "$(date --iso-8601=seconds)" \
        > "$SUCCESS_FLAG"

else

    echo
    echo "ERROR: Required outputs are missing."
    echo "The script can safely be rerun."
    echo "Existing valid outputs will be skipped."

    exit 1

fi


# =========================================================================
# 13. Final report
# =========================================================================

echo
echo "============================================================"
echo "Kraken2 aggregation completed successfully"
echo "============================================================"

echo
echo "Number of samples:"
echo "  ${#KRAKEN_FILES[@]}"

echo
echo "Combined results:"
echo "  $COMBINED_DIR"

echo
echo "Sample summaries:"
echo "  $SUMMARY_DIR"

echo
echo "Abundance tables:"
echo "  $ABUNDANCE_DIR"

echo
echo "Visualisations:"
echo "  $VIS_DIR"

echo
echo "Krona:"
echo "  $KRONA_DIR"

echo
echo "Completion flag:"
echo "  $SUCCESS_FLAG"

echo
echo "============================================================"
echo "Done."
echo "============================================================"

