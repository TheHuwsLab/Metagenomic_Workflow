#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH --job-name=Metaphlan_DM
#SBATCH --error=/users/3057556/jobs/Metaphlan_DM_%j.error
#SBATCH --output=/users/3057556/jobs/Metaphlan_DM_%j.txt
#SBATCH --partition=k2-bioinf
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load Anaconda and activate environment
module load apps/anaconda3/2024.06/bin
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /mnt/scratch2/igfs-anaconda/conda-envs/metaphlan_4.2.5

# Navigate to target directory
cd /mnt/scratch2/users/3057556/DawnMeats/Azenta/40-1311511016/00_fastq/samples/data || exit 1

DIR_TAGS=("K")

for DIR_TAG in "${DIR_TAGS[@]}"; do
  # Look for sample directories matching the tag prefix
  for dir in ./"${DIR_TAG}"*; do
    if [[ -d "$dir" ]]; then

      # Process R1 files located strictly inside the current sample directory
      for file in "$dir"/*_unmapped_trimmed_paired_R1.fastq.gz; do
        # Ensure files actually exist before executing
        [[ -e "$file" ]] || break

        filename=$(basename "$file")
        sample_name=$(echo "$filename" | sed 's/_unmapped_trimmed_paired_R1.fastq.gz//')

        # Build exact paths for forward and reverse reads
        sample_name_1="${dir}/${sample_name}_unmapped_trimmed_paired_R1.fastq.gz"
        sample_name_2="${dir}/${sample_name}_unmapped_trimmed_paired_R2.fastq.gz"

        # Structure the target output directory paths
        outdir="${dir}/${sample_name}_metaphlan"
        bowtie2_output="${outdir}/${sample_name}_metaphlan.bowtie2.bz2"
        report="${outdir}/${sample_name}_metaphlan_report.txt"

        echo "Processing directory: $dir | Sample: $sample_name"
        mkdir -p "$outdir"

        # Run MetaPhlAn 4 targeting both paired fastq files
        metaphlan "${sample_name_1},${sample_name_2}" \
          --input_type fastq \
          --mapout "$bowtie2_output" \
          -o "$report" \
          --nproc 20
      done

    fi
  done
done
