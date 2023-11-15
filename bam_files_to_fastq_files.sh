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

module load SAMtools/1.15-GCC-11.2.0

# Important note: To save storage size, we need to obtain fastq files as fq.gz extension instead of .fq extension.
# So, for this, command below can be used:

samtools fastq -TMM,ML --threads 128 sample_01_merged.bam | gzip > sample_01_merged.fq.gz
exit

# MM and ML tags for methylation are kept in fastq file after conversion of BAM file to fastq file.

# Its also possible to get uncompressed fastq files through command below:

# samtools fastq -TMM,ML --threads 128 sample_01_merged.bam > sample_01_merged.fq

# samtools fastq -T '*' --threads 128 sample_01_merged.bam | gzip > sample_01_merged.fq.gz
# samtools fastq -T '*' --threads 128 sample_01_merged.bam > sample_01_merged.fq

# How to pipe conversion of BAM to FASTQ and align to reference genome: https://github.com/nanoporetech/dorado/issues/145
# If you do: samtools fastq -T '*' you can carry all the tags from the SAM without loss.
# Note: you need at least samtools 1.16 for this piped command below

# samtools fastq -T '*' --threads 128 sample_01_merged.bam | minimap2 -k17 -ax map-ont --secondary=yes -t 128 -y /path_to_directory/UCSC_Mus_musculus.GRCm38_genome.fa.gz - | samtools view -Sb - --threads 128 | samtools sort - --threads 128 > sample_01_merged_aligned.bam
