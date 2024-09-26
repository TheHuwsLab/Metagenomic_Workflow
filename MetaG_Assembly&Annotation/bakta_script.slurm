#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --job-name=&&
#SBATCH --error=/users/&&/jobs/&&.error
#SBATCH --output=/users/&&/jobs/&&.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --mail-user=&&@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL



module load apps/anaconda3/2024.106/bin

conda activate /&&/conda-envs/bakta_1.8.2

# Set location of Bakta database to use:
DBNAME=/&&/conda-envs/bakta_1.8.2/db


cd /&&/data



for dir in ./*; do [[ -d "$dir" ]] && for file in "$dir"/'*_unmapped_metaspades/contigs.fasta'; do sample_name=$(echo $file | sed 's/\.\///;s/\/.*//' ); outdir=$(echo ${sample_name}/${sample_name}_bakta); bakta_output=$(echo ${sample_name}_bakta);  echo "$outdir"; mkdir "$outdir" ; bakta --db $DBNAME --prefix "$bakta_output" --output "$outdir" --min-contig-length 1000  --keep-contig-headers --meta --threads 24 -v --force $file ; done; done