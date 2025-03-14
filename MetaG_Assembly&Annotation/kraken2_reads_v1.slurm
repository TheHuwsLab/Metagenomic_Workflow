#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --job-name=Kraken_Reads_Yasmin_Test
#SBATCH --error=/users/3057556/jobs/Kraken2_Reads_YTest.error
#SBATCH --output=/users/3057556/jobs/Kraken2_Reads_YTest.txt
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL



module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/kraken2_2.1.3

# Set location of Kraken2 database to use:
DBNAME=/mnt/scratch2/igfs-anaconda/conda-dbs/kraken2/k2_pluspfp_20240904


cd /mnt/scratch2/users/3057556/Yasmin/testing


for dir in ./*; do 
    if [[ -d "$dir" ]]; then  
        for file in "$dir"/*_R1_001.fastq.gz; do 
            filename=$(basename "$file")  
            sample_name=$(echo "$file" | sed 's/\.\///;s/\/.*//')  
            sample_name_1="${filename}"  
            sample_name_2=$(echo "$sample_name_1" | sed 's/R1/R2/g')  
            
            outdir="${sample_name}/${sample_name}_kraken2"  
            kraken_output="${outdir}/${sample_name}_kraken2_reads.txt"  
            report="${outdir}/${sample_name}_kraken2_reads_report.txt"  

            echo "Creating output directory: $outdir"
            mkdir -p "$outdir"  

            echo "Running Kraken2 on: $sample_name_1 and $sample_name_2"
            kraken2 --threads 20 \
                    --output "$kraken_output" \
                    --report "$report" \
                    --db "$DBNAME" \
                    --use-names \
                    --paired "$dir/$sample_name_1" "$dir/$sample_name_2"  
        done  
    fi  
done