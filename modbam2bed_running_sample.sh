#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=sssss@sssss.com # for example user@domain.com
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=14-00:00:00 # run for 14 days maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

# Load modules

module load SAMtools/1.15-GCC-11.2.0
module load modbam2bed/0.6.2

# Run modbam2bed

modbam2bed --aggregate -e -m 5mC --cpg -p sample_01_CpGs --threads 128 -c /path_to_directory/UCSC_Mus_musculus.GRCm38_genome.fa.gz sample_01_merged_and_aligned_sorted.bam
exit
