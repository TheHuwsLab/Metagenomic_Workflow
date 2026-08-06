#!/bin/bash

#SBATCH --job-name=Contam_Report_DM
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --error=/users/3057556/jobs/Contam_Report_DM_%j.error
#SBATCH --output=/users/3057556/jobs/Contam_Report_DM_%j.log
#SBATCH --partition=k2-hipri,k2-bioinf
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

JOB_DIR="/mnt/scratch2/users/3057556/DawnMeats/Azenta/40-1311511016/00_fastq/samples/data"
DIR_TAGS=("K")


# Enable nullglob
shopt -s nullglob

# Output Files
OUTPUT_CSV="$JOB_DIR/mapping_summary_test.csv"           # The "Bridge" file (Long IDs)
OVERVIEW_CSV="$JOB_DIR/sample_overview.csv"               # High-level summary
OUTPUT_DB_CSV="$JOB_DIR/mapping_by_db_test_transposed.csv" # The Transposed species table
OUTPUT_DB_TEMP="$JOB_DIR/mapping_by_db_temp.txt"

# Initialize files
echo "Sample_ID,Mapped_Reads,Total_Reads,Percent_Mapped" > "$OUTPUT_CSV"
echo "sample,total_reads,mapped_reads,percent_mapped" > "$OVERVIEW_CSV"
> "$OUTPUT_DB_TEMP"

# Change to working directory
cd "$JOB_DIR" || { echo "Error: Cannot change directory to $JOB_DIR" >&2; exit 1; }

# Process each directory
for DIR_TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"${DIR_TAG}"*; do
    if [[ -d "$dir" ]]; then
      raw_sample_name=$(basename "$dir")
      bam_file=$(ls "$dir"/*_sorted.bam 2>/dev/null | head -n 1 || true)

      if [[ -z "$bam_file" ]]; then
        echo "${raw_sample_name},NA,0,0.00" >> "$OVERVIEW_CSV"
        continue
      fi

      total_reads=$(samtools view -@ 2 -c "$bam_file" 2>/dev/null || echo 0)
      mapped_reads=$(samtools view -@ 2 -c -F 4 "$bam_file" 2>/dev/null || echo 0)

      if [[ "$total_reads" -gt 0 ]]; then
        percent=$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN{printf "%.2f", (m/t)*100}')
      else
        percent="0.00"
      fi

      # 1. Write to Per-Sample Overview
      echo "${raw_sample_name},${total_reads},${mapped_reads},${percent}" >> "$OVERVIEW_CSV"

      # 2. Extract Genome-Specific Stats
      if [[ -f "${bam_file}.bai" ]] || [[ -f "${bam_file%%.bam}.bai" ]]; then
        samtools idxstats "$bam_file" 2>/dev/null | awk -v sample="$raw_sample_name" -v total="$total_reads" -v percent="$percent" '
        {
          if ($1 != "*" && $3 > 0) {
            split($1, a, "_");
            prefix = (length(a) >= 2) ? a[1] "_" a[2] : a[1];
            counts[prefix] += $3;
          }
        }
        END {
          if (length(counts) == 0) {
            print sample, "NO_MAPPING", 0;
          } else {
            for (p in counts) {
              # Format: SampleID, GenomePrefix, Count, TotalReads, PercentMapped
              print sample, p, counts[p], total, percent;
            }
          }
        }' >> "$OUTPUT_DB_TEMP"
      else
        samtools view -F 4 "$bam_file" 2>/dev/null | awk -v sample="$raw_sample_name" -v total="$total_reads" -v percent="$percent" '
        {
          ref = $3;
          if (ref != "*") {
            split(ref, a, "_");
            prefix = (length(a) >= 2) ? a[1] "_" a[2] : a[1];
            counts[prefix]++;
          }
        }
        END {
          if (length(counts) == 0) {
            print sample, "NO_MAPPING", 0, total, percent;
          } else {
            for (p in counts) {
              print sample, p, counts[p], total, percent;
            }
          }
        }' >> "$OUTPUT_DB_TEMP"
      fi
    fi
  done
done

# --- FINAL FILE GENERATION ---

# 1. Generate the Bridge File (mapping_summary_test.csv)
# Creates: "K21-Facility1-A Lolium_perenne 5357,5357,20043552,97.92"
grep -vE "NO_MAPPING" "$OUTPUT_DB_TEMP" | awk '{print $1, $2, $3 "," $4 "," $5}' >> "$OUTPUT_CSV"

# 2. Generate the Transposed Species Table (mapping_by_db_test_transposed.csv)
# This uses a pure AWK pivot to create Sample Rows and Genome Columns.
awk '
{
    if ($2 != "NO_MAPPING") {
        sample=$1; genome=$2; count=$3;
        data[sample][genome] = count;
        samples[sample] = 1;
        genomes[$2] = 1;
    }
}
END {
    # Print Header
    printf "sample";
    # Sort genomes alphabetically for columns
    n = asort(genomes, sorted_genomes);
    for (i=1; i<=n; i++) printf ",%s", sorted_genomes[i];
    print "";

    # Print Rows
    for (s in samples) {
        printf "%s", s;
        for (i=1; i<=n; i++) {
            g = sorted_genomes[i];
            val = (data[s][g] ? data[s][g] : 0);
            printf ",%d", val;
        }
        print "";
    }
}' "$OUTPUT_DB_TEMP" > "$OUTPUT_DB_CSV"

# Clean up
rm -f "$OUTPUT_DB_TEMP"

echo "Processing Complete."
echo "1. High-level summary: $OVERVIEW_CSV"
echo "2. Bridge file (for database join): $OUTPUT_CSV"
echo "3. Transposed species table: $OUTPUT_DB_CSV"