#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --job-name=trimm
#SBATCH --error=Trim_$$.error
#SBATCH --output=Trim_$$.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

conda activate /&&/conda-envs/trimmomatic_0.39/


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
        "$dir/$outfile_1.trimmed.fastq.gz" "$dir/$outfile_1.unpaired.fastq.gz" \
        "$dir/$outfile_2.trimmed.fastq.gz" "$dir/$outfile_2.unpaired.fastq.gz" \
        ILLUMINACLIP:/users/3057556/aux/TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    done
  fi
done
