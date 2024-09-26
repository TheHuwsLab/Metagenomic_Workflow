#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --job-name=$$
#SBATCH --error=/users/$$/jobs/$$.error
#SBATCH --output=/users/$$/jobs/$$.txt
#SBATCH --partition=bio-compute,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=100:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL


module load apps/anaconda3/2024.06/bin

source activate /&&/conda-envs/eggnog-mapper-2.1.12


# Set location of eggnog-mapper database to use:
export EGGNOG_DATA_DIR=/&&/conda-envs/eggnog-mapper-2.1.12/db

cd /mnt/scratch2/users/$$/$$


for dir in ./*; do
    if [[ -d "$dir" ]]; then
        for file in "$dir"/*_pyrodigal/*_pyrodigal_aa.fa; do
            sample_name=$(echo "$file" | sed 's/\.\///;s/\/.*//')
            outdir="${sample_name}/${sample_name}_eggnog_mapper"
            eggnog_output="${sample_name}_pyrodigal_eggnog_mapped"

            # Check if the file already exists
            if [[ -f "${outdir}/${sample_name}_pyrodigal_eggnog_mapped.emapper.annotations.xlsx" ]]; then
                echo "Already computed $sample_name"


            else
                echo "$outdir"
                mkdir -p "$outdir"
                emapper.py -i "$file" --output_dir "$outdir" --output "$eggnog_output" --score 60 --subject_cover 60 --sensmode sensitive --dbmem --decorate_gff yes --excel --cpu 24 --override
            fi
        done
    fi
done
