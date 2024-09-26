#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --job-name=$$
#SBATCH --error=/users/$$/jobs/$$.error
#SBATCH --output=/users/$$/jobs/$$.txt
#SBATCH --partition=bio-compute,k2-medpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL



module load apps/anaconda3/2024.06/bin

conda activate /&&/conda-envs/kraken2_2.1.3

# Set location of Kraken2 database to use:
DBNAME=/&&/conda-dbs/kraken2/2_pluspfp_20230605

cd /mnt/scratch2/users/$$/$$


for dir in ./*; do
  if [[ -d "$dir" ]]; then
    for file in "$dir"/*_unmapped_metaspades/contigs.fasta; do
      sample_name=$(echo "$file" | sed 's/\.\///;s/\/.*//')
      outdir="${sample_name}/${sample_name}_kraken2"
      kraken_output="${outdir}/${sample_name}_kraken2.txt"
      report_mpa="${outdir}/${sample_name}_kraken2_report_mpa.txt"
      report="${outdir}/${sample_name}_kraken2_report.txt"

      echo "$outdir"
      mkdir -p "$outdir"
      
      kraken2 --threads 24 --output "$kraken_output" --report "$report_mpa" \
        --db "$DBNAME" --use-names --use-mpa-style "$file"
      
      kraken2 --threads 24 --output "$kraken_output" --report "$report" \
        --db "$DBNAME" --use-names "$file"
    done
  fi
done
