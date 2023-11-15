#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH mail-user=sssss@sssss.com # for example user@domain.com
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=1 # number of nodes [mimir 66 ve 72 arasında bir node atayacak rastgele]
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=14-00:00:00 # run for 4 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

ml load SAMtools/1.15-GCC-11.2.0

# Increasing how many file can be read at the same time

# -Sn 50000 or

ulimit -n 5000

# samtools merge to merge bam files into one larger bam file

samtools merge -n -c -@ 128 -o sample_01_merged.bam bam_pass/*.bam
exit

# interactive run: $ srun --job-name "sample01" --cpus-per-task 1 --mem-per-cpu 3900 bash samtools_merge_BAM_files.sh
