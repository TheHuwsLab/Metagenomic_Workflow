#!/bin/bash

#SBATCH --cpus-per-task=30
#SBATCH --mem=120G
#SBATCH --job-name=Yasmin_contam_test
#SBATCH --error=/users/3057556/jobs/contam_yasmin_test.error
#SBATCH --output=/users/3057556/jobs/contam_yasmin_test.txt
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/bowtie2_2.5.1

cd /mnt/scratch2/users/3057556/Yasmin/testing


for dir in ./*; do
  if [[ -d "$dir" ]]; then
    for file in "$dir"/*R1_001.fastq.gz; do
      filename=$(basename "$file")
      sample_name_1="${filename}"
      sample_name_2=$(echo "$sample_name_1" | sed 's/R1/R2/g')

      # Define output filenames
      unmapped_sorted_bam=$(echo "$sample_name_1" | sed 's/_R1_001.fastq.gz/_sorted_unmapped.bam/g')
      unmapped_R1=$(echo "$sample_name_1" | sed 's/_R1_001.fastq.gz/_mapped_R1.fastq/g')
      unmapped_R2=$(echo "$sample_name_1" | sed 's/_R1_001.fastq.gz/_mapped_R2.fastq/g')

      # Run bowtie2 and pipe the output directly into samtools without writing intermediate files
      bowtie2 --threads 30 --very-sensitive-local \
        -x /mnt/scratch2/igfs-databases/reference_genomes/combined_cow_sheep_human/Human_Sheep_Cow_Combined_Bowtie2 \
        -1 "$dir/$sample_name_1" -2 "$dir/$sample_name_2" | \
      samtools view -b -f 12 - | \
      samtools sort -n - -o "$dir/$unmapped_sorted_bam"

      # Convert unmapped reads to FASTQ
      bedtools bamtofastq -i "$dir/$unmapped_sorted_bam" -fq "$dir/$unmapped_R1" -fq2 "$dir/$unmapped_R2"

    done
  fi
done
