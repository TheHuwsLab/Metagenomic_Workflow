#!/bin/bash

#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --job-name=Py_$$
#SBATCH --error=/users/&&/jobs/Pyrodigal_$$.error
#SBATCH --output=/users/&&/jobs/Pyrodigal_$$.txt
#SBATCH --partition=&&
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --mail-user=$$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin
source activate /&&/conda-envs/pyrodigal_3.0.1

cd /&&/data

# Define the directory tag for sample directories
DIR_TAG="PN"

# Loop through all directories in the current folder
for dir in "${DIR_TAG}"*; do
  # Check if the element is a directory
  if [[ -d "$dir" ]]; then
    # Loop through all files matching the pattern in the subdirectories
    for file in "$dir"/*_unmapped_metaspades/contigs.fasta; do
      # Extract the sample name from the file path
      sample_name=$(echo $file | sed 's/\.\///;s/\/.*//')
      # Define the output directory and output file paths
      outdir="${sample_name}/${sample_name}_pyrodigal"
      pyrodigal_out_gff="${outdir}/${sample_name}_pyrodigal.gff"
      pyrodigal_out_nt="${outdir}/${sample_name}_pyrodigal_nt.fa"
      pyrodigal_out_aa="${outdir}/${sample_name}_pyrodigal_aa.fa"

      # Check if all output files exist and are not empty
      if [[ -s "$pyrodigal_out_gff" && -s "$pyrodigal_out_nt" && -s "$pyrodigal_out_aa" ]]; then
        echo "All output files for $sample_name exist and are not empty. Skipping."
        continue
      fi

      echo "$outdir"
      # Create the output directory
      mkdir "$outdir"
      # Run Pyrodigal with the specified input and output files
      pyrodigal -i "$file" \
                -o "$pyrodigal_out_gff" \
                -a "$pyrodigal_out_aa" \
                -d "$pyrodigal_out_nt" \
                -p meta -j 12
    done
  fi
done
