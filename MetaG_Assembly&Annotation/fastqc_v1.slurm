#!/bin/bash

#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --job-name=fastqc
#SBATCH --error=/users/&&/jobs/&&.error
#SBATCH --output=/users/&&/jobs/&&.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mail-user=&&
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

conda activate /&&/conda-envs/multiqc

# The working directory to where your samples are
search_dir=/&&/data

# Set the path to the directory where you want to store FastQC results
output_dir=/&&/FastQC_Filtered

mkdir $output_dir

# Find all FASTQ files in subdirectories and loop through them

find $search_dir -type f -name "*unmapped_R*.fastq" | xargs fastqc -t 12 -o "$output_dir"
