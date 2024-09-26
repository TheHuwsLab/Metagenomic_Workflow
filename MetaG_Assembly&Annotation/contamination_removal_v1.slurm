#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --job-name=$$
#SBATCH --error=/users/$$/jobs/$$.error
#SBATCH --output=/users/$$/jobs/$$.txt
#SBATCH --partition=bio-compute
#SBATCH --nodes=1
#SBATCH --time=200:00:00
#SBATCH --mail-user=$$@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

conda activate /&&/conda-envs/bowtie2_2.5.1

cd /mnt/scratch2/users/$$/$$


for dir in ./*; do
  if [[ -d "$dir" ]]; then
    for file in "$dir"/*R1_001.trimmed.fastq.gz; do
      filename=$(basename "$file")
      sample_name_1="${filename}"
      sample_name_2=$(echo "$sample_name_1" | sed 's/R1/R2/g')
      mapped_bam=$(echo "$sample_name_1" | sed 's/_R1_001.trimmed.fastq.gz/_mapped.bam/g')
      unmapped_bam=$(echo "$mapped_bam" | sed 's/mapped/unmapped/g')
      unmapped_sorted_bam=$(echo "$unmapped_bam" | sed 's/unmapped/sorted_unmapped/g')
      unmapped_R1=$(echo "$unmapped_bam" | sed 's/mapped/mapped_R1/g' | sed 's/.bam/.fastq/g')
      unmapped_R2=$(echo "$unmapped_bam" | sed 's/mapped/mapped_R2/g' | sed 's/.bam/.fastq/g')

      bowtie2 --threads 20 -x /mnt/scratch/igfs-databases/reference_genomes/combined_cow_sheep_human/Human_Sheep_Cow_Combined_Bowtie2/Human_Sheep_Cow_Combined \
        -1 "$dir/$sample_name_1" -2 "$dir/$sample_name_2" | samtools view -Sbh - > "$dir/$mapped_bam"
      
      samtools view -b -f 12 "$dir/$mapped_bam" > "$dir/$unmapped_bam"
      samtools sort -n "$dir/$unmapped_bam" -o "$dir/$unmapped_sorted_bam"
      bedtools bamtofastq -i "$dir/$unmapped_sorted_bam" -fq "$dir/$unmapped_R1" -fq2 "$dir/$unmapped_R2"
    done
  fi
done
