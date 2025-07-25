#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=build_rumen_bowtie2
#SBATCH --error=/users/$$/jobs/build_rumen_bowtie2_error.log
#SBATCH --output=/users/$$/jobs/build_rumen_bowtie2_output.log
#SBATCH --partition=$$
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load the required modules
module load apps/anaconda3/2024.06/bin
source activate /mnt/$$/conda-envs/bowtie2_2.5.1

cd /mnt/$$/reference_genomes/rumen_contamination_genomes

bowtie2-build rumen_contamination_genomes.fasta  rumen_contamination_bt2_db --threads 20 --large-index
