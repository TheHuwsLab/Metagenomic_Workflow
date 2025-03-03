#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=ReadCov_Sophie_CP001719_1
#SBATCH --error=/users/3057556/jobs/read_coverage/Sophie_CP001719_1_error.log
#SBATCH --output=/users/3057556/jobs/read_coverage/Sophie_CP001719_1_output.log
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

WORK_DIR="/mnt/scratch2/users/3057556/SOPHIE/FastQ/CP001719_1_Readmapped"

# Change to the working directory
cd "$WORK_DIR" || exit

# Define the output file
output_file="CP001719_1_read_coverage_stats.csv"

# Write header to the CSV file
echo -e "Sample\tTotal_Reads\tMapped_Reads\tMean_Coverage\tNormalised_Coverage" > "$output_file"

# Loop through each BAM file
for bam in *.bam; do
    if [[ -f "$bam" ]]; then
        sample_name=$(basename "$bam" .bam)

        # Get total reads in BAM (including unmapped)
        total_reads=$(samtools view -c "$bam")

        # Get total mapped reads (excluding unmapped reads)
        mapped_reads=$(samtools view -c -F 4 "$bam")

        # Get coverage statistics using bedtools
        mean_coverage=$(bedtools genomecov -ibam "$bam" -d | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')

        if [[ "$total_reads" -gt 0 ]]; then
            # Calculate normalised coverage
            normalised_coverage=$(echo "$mean_coverage / $total_reads" | bc -l)

            # Check if the value has more than 6 decimal places
            if (( $(echo "$normalised_coverage < 0.000001" | bc -l) )); then
                # Use scientific notation if < 0.000001
                normalised_coverage=$(printf "%.3e" "$normalised_coverage")
            else
                # Otherwise, round to 6 decimal places
                normalised_coverage=$(printf "%.6f" "$normalised_coverage")
            fi
        else
            normalised_coverage="0.000000"
        fi

        # Append results to the output file
        echo -e "${sample_name}\t${total_reads}\t${mapped_reads}\t${mean_coverage}\t${normalised_coverage}" >> "$output_file"
    fi
done

echo "Report saved as $output_file"
