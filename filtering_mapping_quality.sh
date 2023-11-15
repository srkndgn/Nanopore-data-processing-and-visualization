#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=sssss@sssss.com # for example user@domain.com
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=14-00:00:00 # run for 6 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

# Load samtools module on Elja

module load SAMtools/1.15-GCC-11.2.0

# Filter reads with less than mapping quality 30

# Define an array of BAM file names
bam_files=(
    "/path/to/bam/file1.bam"
    "/path/to/bam/file2.bam"
    "/path/to/bam/file3.bam"
)

# Loop through the list of BAM files
for bam_file in "${bam_files[@]}"
do
    samtools view -h -b -q 30 --threads 128 "$bam_file" > "${bam_file%.*}_mapping_quality_30_filtered.bam"
done
