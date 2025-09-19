#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --job-name=trimm
#SBATCH --error=/users/&&/jobs/trimm_&&.error
#SBATCH --output=/users/&&/jobs/trimm_&&.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL


module load apps/anaconda3/2024.06/bin

# Activate the conda environment with Trimmomatic installed
source activate /&&/conda-envs/trimmomatic_0.39

# Change the working directory to where your samples are
cd /&&/data

for dir in ./*; do
  if [[ -d "$dir" ]]; then
    for file in "$dir"/*R1_001.fastq.gz; do
      filename=$(basename "$file")
      sample_name_1="${filename}"
      sample_name_2=$(echo "$sample_name_1" | sed 's/R1/R2/g')
      outfile_1=$(echo "$sample_name_1" | sed 's/.fastq.gz//g')
      outfile_2=$(echo "$sample_name_2" | sed 's/.fastq.gz//g')

      trimmomatic PE -phred33 -threads 8 \
        "$dir/$sample_name_1" "$dir/$sample_name_2" \
        "$dir/$outfile_1_trimmed_paired.fastq.gz" "$dir/$outfile_1_trimmed_unpaired.fastq.gz" \
        "$dir/$outfile_2_trimmed_paired.fastq.gz" "$dir/$outfile_2_trimmed_unpaired.fastq.gz" \
        ILLUMINACLIP:$$/conda-envs/trimmomatic_0.39/adaptors/TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    done
  fi
done
