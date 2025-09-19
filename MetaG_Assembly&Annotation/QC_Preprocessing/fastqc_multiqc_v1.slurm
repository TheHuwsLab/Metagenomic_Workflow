#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --job-name=fastqc_mqc_$$
#SBATCH --error=/users/&&/jobs/fastqc_mqc_&&.error
#SBATCH --output=/users/&&/jobs/fastqc_mqc_&&.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=&&
#SBATCH --mail-type=BEGIN,END,FAIL

# Activate the conda environment with FastQC installed
module load apps/anaconda3/2024.06/bin
source activate /&&/conda-envs/multiqc_1.30

# The working directory to where your samples are
search_dir=/&&/data

# Set the path to the directory where you want to store FastQC results
fqc_output_dir=/&&/FastQC
mqc_output_dir=/&&/MultiQC

mkdir -p $fqc_output_dir
mkdir -p $mqc_output_dir

# Find all FASTQ files in subdirectories and loop through them
find $search_dir -type f -name "*R*_001.fastq.gz" | xargs fastqc -t 8 -o "$fqc_output_dir"
multiqc --outdir $mqc_output_dir $fqc_output_dir




