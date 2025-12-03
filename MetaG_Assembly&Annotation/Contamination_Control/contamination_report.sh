#!/bin/bash

#SBATCH --job-name=Contam_Report_GL3
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --error=/users/3057556/jobs/Contam_Report_GL3_%j.error
#SBATCH --output=/users/3057556/jobs/Contam_Report_GL3_%j.log
#SBATCH --partition=k2-hipri
#SBATCH --time=03:00:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

JOB_DIR="/mnt/scratch2/users/3057556/data/"
DIR_TAGS=("E" "L" "P")

# Enable nullglob to handle case where no directories match pattern
shopt -s nullglob

OUTPUT_CSV="$JOB_DIR/mapping_summary_test.csv"
OUTPUT_DB_CSV="$JOB_DIR/mapping_by_db_test.csv"
OUTPUT_DB_TEMP="$JOB_DIR/mapping_by_db_temp.txt"

# Initialize output files with headers
echo "sample,sorted_bam,total_reads,mapped_reads,percent_mapped" > "$OUTPUT_CSV"

# Clear temp file for collecting genome data
> "$OUTPUT_DB_TEMP"

# Change to working directory
cd "$JOB_DIR" || { echo "Error: Cannot change directory to $JOB_DIR" >&2; exit 1; }

# Process each directory matching the pattern for each tag
for DIR_TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"${DIR_TAG}"*; do
    if [[ -d "$dir" ]]; then
      sample_name=$(basename "$dir")

      # Find the sorted BAM file
      bam_file=$(ls "$dir"/*_sorted.bam 2>/dev/null | head -n 1 || true)

      # Handle missing BAM file
      if [[ -z "$bam_file" ]]; then
        echo "${sample_name},NA,0,0,0.00" >> "$OUTPUT_CSV"
        echo "${sample_name} NO_BAM 0" >> "$OUTPUT_DB_TEMP"
        continue
      fi

      # Calculate total and mapped reads
      total_reads=$(samtools view -@ 2 -c "$bam_file" 2>/dev/null)
      if [[ $? -ne 0 ]]; then
        total_reads=0
      fi

      mapped_reads=$(samtools view -@ 2 -c -F 4 "$bam_file" 2>/dev/null)
      if [[ $? -ne 0 ]]; then
        mapped_reads=0
      fi

      # Calculate percentage
      if [[ "$total_reads" -gt 0 ]]; then
        percent=$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN{printf "%.2f", (m/t)*100}')
      else
        percent="0.00"
      fi

      # Write summary statistics
      echo "${sample_name},$(basename "$bam_file"),$total_reads,$mapped_reads,$percent" >> "$OUTPUT_CSV"

      # Aggregate mapped reads by genome (database sequence prefix up to second '_')
      # Check if BAM index exists (more efficient than running idxstats to test)
      if [[ -f "${bam_file}.bai" ]] || [[ -f "${bam_file%%.bam}.bai" ]]; then
        # Use idxstats (fast method with index)
        samtools idxstats "$bam_file" 2>/dev/null | awk -v sample="$sample_name" '
        {
          # Skip unmapped reads (marked with "*" as reference name)
          if ($1 != "*" && $3 > 0) {
            # Extract genome prefix (up to second underscore)
            split($1, a, "_");
            if (length(a) >= 2) {
              prefix = a[1] "_" a[2];
            } else {
              prefix = a[1];
            }
            counts[prefix] += $3;
          }
        }
        END {
          # Output results in space-delimited format for later processing
          if (length(counts) == 0) {
            print sample, "NO_MAPPING", 0;
          } else {
            for (p in counts) {
              print sample, p, counts[p];
            }
          }
        }' >> "$OUTPUT_DB_TEMP"
      else
        # Fallback: Parse mapped reads directly (slower, no index required)
        samtools view -F 4 "$bam_file" 2>/dev/null | awk -v sample="$sample_name" '
        {
          ref = $3;
          # Skip unmapped reads (should not occur with -F 4, but be safe)
          if (ref != "*") {
            # Extract genome prefix (up to second underscore)
            split(ref, a, "_");
            if (length(a) >= 2) {
              prefix = a[1] "_" a[2];
            } else {
              prefix = a[1];
            }
            counts[prefix]++;
          }
        }
        END {
          # Output results in space-delimited format for later processing
          if (length(counts) == 0) {
            print sample, "NO_MAPPING", 0;
          } else {
            for (p in counts) {
              print sample, p, counts[p];
            }
          }
        }' >> "$OUTPUT_DB_TEMP"
      fi
    fi
  done
done

# Clean up temp file
rm -f "$OUTPUT_DB_TEMP"

echo "Mapping summary written to $OUTPUT_CSV"
echo "Per-genome summary (transposed) written to $OUTPUT_DB_CSV"