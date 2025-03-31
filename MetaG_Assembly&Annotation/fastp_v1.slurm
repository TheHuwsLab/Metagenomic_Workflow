#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --job-name=fastp_test
#SBATCH --error=/users/3057556/jobs/fastp_test.error
#SBATCH --output=/users/3057556/jobs/fastp_test.txt
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/fastp_0.24.0


#!/bin/bash

#!/bin/bash

# Define the parent directory containing subdirectories
INPUT_DIR="/mnt/scratch2/users/3057556/testing/"

# Loop through each subdirectory in the input directory
for subdir in "$INPUT_DIR"/*/; do
    if [[ -d "$subdir" ]]; then
        echo "Processing directory: $subdir"

        # Detect R1 and R2 files - May need to modify the name catches here
        R1=$(find "$subdir" -type f -name "*R1_001.fastq.gz" | head -n 1)
        R2=$(find "$subdir" -type f -name "*R2_001.fastq.gz" | head -n 1)

        if [[ -f "$R1" && -f "$R2" ]]; then
            # Define output filenames in the same directory
            out_R1="${R1/.fastq.gz/_trimmed_paired.fastq.gz}"
            out_R2="${R2/.fastq.gz/_trimmed_paired.fastq.gz}"
            out_unpaired_R1="${R1/.fastq.gz/_trimmed_unpaired.fastq.gz}"
            out_unpaired_R2="${R2/.fastq.gz/_trimmed_unpaired.fastq.gz}"
            report_html="${subdir}fastp_report.html"

            echo "Read files to be trimmed: \"$R1\" & \"$R2\""

            # Run fastp with unpaired output enabled
            fastp -i "$R1" -I "$R2" -o "$out_R1" -O "$out_R2" \
                  --unpaired1 "$out_unpaired_R1" --unpaired2 "$out_unpaired_R2" \
                  --html "$report_html"  --detect_adapter_for_pe \
                  --thread 20 --qualified_quality_phred 20 --length_required 50

            echo "Finished processing: $(basename "$subdir")"
        else
            echo "WARNING: Missing R1 or R2 files in $subdir, skipping..."
        fi
    fi
done

echo "All samples processed!"





