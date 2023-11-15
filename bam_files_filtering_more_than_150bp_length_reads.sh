#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssssss@ssssss.com # for example user@domain.com
#SBATCH --partition=mimir  # request node from a specific partition
#SBATCH --nodes=1                 # number of nodes
#SBATCH --ntasks-per-node=20      # 48 cores per node (96 in total)
#SBATCH --mem-per-cpu=3900        # MB RAM per cpu core
#SBATCH --time=14-00:00:00         # run for 4 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

ml load SAMtools/1.15-GCC-11.2.0

# Define an array of BAM file names
bam_files=(
    "/path/to/bam/file1.bam"
    "/path/to/bam/file2.bam"
    "/path/to/bam/file3.bam"
)

# Loop through the list of BAM files
for bam_file in "${bam_files[@]}"
do
    samtools view -@ 128 -h "$bam_file" | awk 'length($10) >= 150 || $1 ~ /^@/' | samtools view -@ 128 -b -o "${bam_file%.*}_150bp_read_length_filtered.bam" -
done
