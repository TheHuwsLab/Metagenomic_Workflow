#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=multiqc
#SBATCH --error=MultiQC_$$.error
#SBATCH --output=MultiQC_$$.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

conda activate /&&/conda-envs/multiqc

# Change the working directory to where your FastQC output files are
cd /&&/data
# Set the path to the directory where you want to store FastQC results
output_dir=/path/to/output_directory

# Find all FASTQ files in subdirectories and loop through them
multiqc --outdir $output_dir .
