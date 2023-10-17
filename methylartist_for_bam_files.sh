
# Methylartist: tools for visualizing modified bases from nanopore sequence data
# Tools for parsing and plotting methylation patterns from nanopore data
# Reference page > https://github.com/adamewing/methylartist 
# Citation > Seth W Cheetham, Michaela Kindlova, Adam D Ewing, Methylartist: tools for visualizing modified bases from nanopore sequence data, Bioinformatics, Volume 38, Issue 11, 1 June 2022, Pages 3109–3112, https://doi.org/10.1093/bioinformatics/btac292


# To execute the contents of the .bashrc file in the current shell session (necessary to start conda)
    source ~/.bashrc

# In a clean environment:

    cd anaconda3/envs/                                   # go to your anaconda3 directory
    
    conda create -n methylartist python=3.8 anaconda     # create methylartist environment and python
    
    cd anaconda3/envs/methylartist/                       # go to your working environment
    
    conda activate methylartist                           # activate conda environment
    
    pip install methylartist                              # methylartist

# create your working directory for your methylartist environment
    mkdir methylartist             # working directory should be outside the anaconda3 directory
    cd methylartist                # go to your modbamtools working directory

# activate the methylartist environment in your working directory
    conda activate methylartist             #activate conda environment

# Download chr.gtf file for mm10
https://ftp.ensembl.org/pub/release-91/gtf/mus_musculus/


# To sort chr.gtf files
# To use tabix (Generic indexer for TAB-delimited genome position files) install Samptools, bcftools and HTSlib 
# use shark window to load thse modules

module load library/htslib/1.16/gcc-8.5.0
module load genomics/ngs/bcftools/1.11/gcc-8.3.1
module load genomics/ngs/samtools/1.16.1/gcc-8.5.0

# To check installed modules
module li

#############################################################
    # gtf.gz to .gtf.gz.tbi
    # convert a GTF file (a gene annotation file) to a tabix index file (a compressed and indexed file for fast retrieval of genomic regions)
    # If your input file is not sorted by chromosome and position, which is required by tabix1. You can check if your file is sorted by using the sort command with the -c option:
#############################################################

zcat Mus_musculus.GRCm38.91.chr.gtf.gz | sort -c -k1,1 -k4,4n -k5,5n -t$'\t'
    
    # If the output says sort: disorder: …, then you need to sort your file by chromosome and position before compressing and indexing it with tabix2. You can do this by using the following commands:

zcat mm10.refGene.gtf.gz | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip > mm10.refGene.sorted.gtf.gz

tabix -p gff mm10.refGene.sorted.gtf.gz

samtools faidx mm10.refGene.sorted.gtf.gz

###############################################################
# To sort bam files
samtools sort sample_01.bam -o sample_01.sorted.bam
samtools sort sample_02.bam -o sample_02.sorted.bam


# To create index file from a sorted bam file
samtools index sample_01.sorted.bam  sample_01.sorted.bam.bai
samtools index sample_02.sorted.bam  sample_02.sorted.bam.bai

###############################################################

# Compressing Files with gzip
# To compress a single file invoke the gzip command followed by the filename:
    gzip filename

# Unzipping gz File
# The command will restore the compressed file to its original state and remove the .gz file
    gzip -d file.gz

# To keep the compressed file pass the -k option to the command:
    gzip -dk file.gz

# To open a .gz file with gunzip simply pass the file name to the command
    gunzip file.gz

# To open and keep both files > .txt and a .gz files with gunzip simply pass the file name to the command
    gunzip file.gz    

# To see the first 10-15 lines of a file in Linux, you can use the head command with the -n option to specify the number of lines to display. Here's the command you can use:
    head -n 15 gencode.vM22.annotation.sorted.gtf
    less gencode.vM22.annotation.sorted.gtf

# To open the file "gencode.vM22.annotation.sorted.gtf" in the nano text editor and view its contents, you can use the following command:
    nano gencode.vM22.annotation.sorted.gtf

# Use the grep command to search for specific patterns or lines in the file:
    grep "pattern" gencode.vM22.annotation.sorted.gtf
# To quit from a file q/ctrl+c/ctrl+x

################################################################

# locus
# Generates smoothed methylation profiles across specific loci with many configurable parameters for one or more samples.

# Gene clusters

methylartist locus -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -i chr18:36771904-37909338 -r GRCm38.p6.genome.fa -g gencode.vM22.annotation.sorted.gtf.gz -n CG  -p 3,5,2,3,3 --width 30 --height 16 --svg
methylartist locus -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -i chr18:36771904-37909338 -r GRCm38.p6.genome.fa -g gencode.vM22.annotation.sorted.gtf.gz -n CG  -p 3,5,2,3,3 --width 30 --height 16

################################################################

# region
# More tractable than "locus" for larger regions using binned methylation.

# Gene clusters

# all with genes names
methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --labelgenes --skip_align_plot --panelratio 7,0,3,7 --svg
methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --labelgenes --skip_align_plot --panelratio 7,0,3,7

# all without genes names
methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --skip_align_plot --panelratio 7,0,3,7 --svg
methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --skip_align_plot --panelratio 7,0,3,7

# just genes with names
methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --labelgenes --skip_align_plot --panelratio 10,0,1,0 --svg
methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --labelgenes --skip_align_plot --panelratio 10,0,1,0

# just genes without names
methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --skip_align_plot --panelratio 10,0,1,0 --svg
methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --skip_align_plot --panelratio 10,0,1,0

################################################################

# You can run all commands with sbatch
sbatch methylartist_all_commands_for_bam.sh

