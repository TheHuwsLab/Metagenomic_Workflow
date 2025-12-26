#!/bin/bash

#SBATCH --job-name=Prepare_E-Mapper_Jobs
#SBATCH --error=/users/3057556/jobs/Prepare_E-Mapper_Jobs_GL4_%j.error
#SBATCH --output=/users/3057556/jobs/Prepare_E-Mapper_Jobs_GL4_%j.txt
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00  # Short time for job submission script
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
module load apps/anaconda3/2024.06/bin

# Define paths
WORK_DIR="/mnt/scratch2/users/3057556/data/"
JOB_DIR="/users/3057556/jobs/GL4_emapper_jobs"


# Create a directory for job scripts if it doesn't exist
mkdir -p "$JOB_DIR"
# Change to the working directory
cd "$WORK_DIR" || exit

# Define the directory tag for sample directories
DIR_TAGS=("E" "L" "P")

# Loop through each sample directory
for TAG in "${DIR_TAGS[@]}"; do
  for dir in ./"$TAG"*; do
    if [[ -d "$dir" ]]; then
      for file in "$dir"/*_pyrodigal/*_pyrodigal_aa.fa; do
        sample_name=$(basename "$dir")
        outdir="${dir}/${sample_name}_eggnog_mapper"
        eggnog_output="${sample_name}_pyrodigal_eggnog_mapped"


        # Check if output files exist and are non-empty
        if [[ -s "${outdir}/${sample_name}_pyrodigal_eggnog_mapped.emapper.annotations.xlsx" ]]; then
          echo "Skipping $sample_name - eggnogMapper results already exist."
          continue
        fi

        # Create job script
        job_script="$JOB_DIR/${sample_name}_eggnogMapper_job.sh"

        cat <<EOL > "$job_script"
#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --job-name=EggnogMapper_${sample_name}
#SBATCH --error=$JOB_DIR/${sample_name}_error.log
#SBATCH --output=$JOB_DIR/${sample_name}_output.log
#SBATCH --partition=k2-medpri,k2-bioinf,k2-lowpri
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=n.dimonaco@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/igfs-anaconda/conda-envs/eggnog-mapper-2.1.12

# Set location of eggnog-mapper database to use:
export EGGNOG_DATA_DIR=/mnt/scratch2/igfs-anaconda/conda-envs/eggnog-mapper-2.1.12/db

mkdir -p "$outdir"

emapper.py -i "$file" --output_dir "$outdir" --output "$eggnog_output" --score 60 \
  --subject_cover 60 --sensmode sensitive --dbmem --decorate_gff yes --excel --cpu 20 --override

# create a success flag only if emapper succeeded
touch "${outdir}/${sample_name}_eggnogmapper.success"

EOL

      # Submit the job
      echo "Submitting job for $sample_name..."
      sbatch "$job_script"

done
    fi
  done
done
