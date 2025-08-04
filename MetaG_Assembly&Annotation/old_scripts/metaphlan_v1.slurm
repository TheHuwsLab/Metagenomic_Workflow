#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --job-name=&&
#SBATCH --error=/users/&&/jobs/Metaphlan_&&.error
#SBATCH --output=/users/&&/jobs/Metaphlan_&&.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --mail-user=&&
#SBATCH --mail-type=BEGIN,END,FAIL



module load apps/anaconda3/2024.06/bin

conda activate /&&/conda-envs/metaphlan_4.0.6

cd /&&/data


for file in *_R1_001.fastq.gz; do
  filename=$(basename "$file")
  sample_name=$(echo "$filename" | sed 's/_R1_001.fastq.gz//')
  
  sample_name_1="${sample_name}_R1_001.fastq.gz"
  sample_name_2="${sample_name}_R2_001.fastq.gz"
  
  outdir="${sample_name}/${sample_name}_metaphlan"
  bowtie2_output="${outdir}/${sample_name}_metaphlan.bowtie2.bz2"
  report="${outdir}/${sample_name}_metaphlan_report.txt"
  
  echo "$outdir"
  mkdir -p "$outdir"
  
  metaphlan "$sample_name_1","$sample_name_2" --bowtie2out "$bowtie2_output" --input_type fastq -o "$report" --nproc 24
done

