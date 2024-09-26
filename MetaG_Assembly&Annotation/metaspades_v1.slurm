#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --job-name=&&
#SBATCH --error=/users/&&/jobs/Metaphlan_&&.error
#SBATCH --output=/users/&&/jobs/Metaphlan_&&.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=200:00:00
#SBATCH --mail-user=&&
#SBATCH --mail-type=BEGIN,END,FAIL



module load apps/anaconda3/2024.06/bin

conda activate /&&/conda-envs/metaphlan_4.0.6

cd /&&/data


for dir in ./*; do
  if [[ -d "$dir" ]]; then
    for file in "$dir"/*_unmapped_R1.fastq; do
      filename=$(basename "$file")
      sample_name=$(echo "$file" | sed 's/\.\///;s/\/.*//')
      sample_name_1="${filename}"
      sample_name_2=$(echo "$sample_name_1" | sed 's/R1/R2/g')
      
      outdir="${sample_name}/${sample_name}_metaphlan"
      bowtie2_output="${sample_name}/${sample_name}_metaphlan/${sample_name}_metaphlan.bowtie2.bz2"
      report="${sample_name}/${sample_name}_metaphlan/${sample_name}_metaphlan_report.txt"
      
      echo "$outdir"
      mkdir -p "$outdir"
      metaphlan "$dir/$sample_name_1","$dir/$sample_name_2" --bowtie2out "$bowtie2_output" --input_type fastq -o "$report" --nproc 24
    done
  fi
done
