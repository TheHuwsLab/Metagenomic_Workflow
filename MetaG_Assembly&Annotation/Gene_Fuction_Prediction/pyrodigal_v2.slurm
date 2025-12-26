#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --job-name=Pyrodigal_GL3
#SBATCH --error=/users/3057556/jobs/Pyrodigal_GL3_%j.error
#SBATCH --output=/users/3057556/jobs/Pyrodigal_GL3_j%.txt
#SBATCH --partition=k2-bioinf,k2-medpri,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/pyrodigal_3.6.3

cd /mnt/scratch2/users/3057556/GLESNI/X204SC24084950-Z01-F003/01.RawData

# Define the directory tag for sample directories
DIR_TAGS=("E" "L" "P")

# Loop through all directories in the current folder
for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    # Check if the element is a directory
    if [[ -d "$dir" ]]; then
      # Loop through all files matching the pattern in the subdirectories
      for file in "$dir"/*_metaspades/contigs.fasta; do
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
                  -p meta -j 10
done
    fi
  done
done

