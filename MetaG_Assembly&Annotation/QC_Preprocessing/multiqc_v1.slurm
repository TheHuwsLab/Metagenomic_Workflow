#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=multiqc
#SBATCH --error=/users/&&/jobs/multi_fastqc_&&.error
#SBATCH --output=/users/&&/jobs/multi_fastqc_&&.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Activate the conda environment with FastQC installed
module load apps/anaconda3/2024.06/bin
source activate /&&/conda-envs/multiqc_1.30

# Change the working directory to where your FastQC output files are
cd /&&/FastQC
# Set the path to the directory where you want to store FastQC results
output_dir=/path/to/output_directory

# Find all FASTQ files in subdirectories and loop through them
multiqc --outdir $output_dir .
