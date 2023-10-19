# modbamtools is a set of tools to manipulate and visualize DNA/RNA base modification data that are stored in bam format
# htslib has included a support for parsing modified base tags from alignment files (MM and ML). 
# These tags have provided a better/efficient way for storing modification data inside alignment files. 
# Reference page > https://rrazaghi.github.io/modbamtools/
# How to cite modbamtools > Roham Razaghi, Paul W. Hook, Shujun Ou, Michael C. Schatz, Kasper D. Hansen, Miten Jain, Winston Timp, Modbamtools: Analysis of single-molecule epigenetic data for long-range profiling, heterogeneity, and clustering
# https://www.biorxiv.org/content/10.1101/2022.07.07.499188v1


# To execute the contents of the .bashrc file in the current shell session (necessary to start conda)
    source ~/.bashrc

# How to view contents of tar file without extracting it
    tar -ztvf my-data.tar.gz
    tar -tvf my-data.tar.gz
    tar -tvf my-data.tar.gz 'search-pattern'


# To extract a tar.gz file
    tar -xf archive.tar.gz

# Install Python
# Required: Python >=3.8

In a clean environment:

    cd anaconda3/envs/                                   # go to your anaconda3 directory
    
    conda create -n modbamtools python=3.8 anaconda      # create Modbamtools environment and install your Python version
    
    cd anaconda3/envs/modbamtools/                       # go to your working environment
    
    conda activate modbamtools                           # activate conda environment
    
    pip install modbamtools                              # modbamtools

# create your working directory for your modbamtools environment
    mkdir modbamtools             # working directory should be outside the anaconda3 directory
    cd modbamtools                # go to your modbamtools working directory

# activate the modbamtools environment in your working directory
    conda activate modbamtools             #activate conda environment


# To sort annotation files
# To use tabix (Generic indexer for TAB-delimited genome position files) install Samptools, bcftools and HTSlib 
# use shark window to load thse modules

module add module load library/htslib/1.16/gcc-8.5.0
module load genomics/ngs/bcftools/1.11/gcc-8.3.1
module load genomics/ngs/samtools/1.16.1/gcc-8.5.0

# To check installed modules
module li

#############################################################
    # gtf.gz to .gtf.gz.tbi
    # convert a GTF file (a gene annotation file) to a tabix index file (a compressed and indexed file for fast retrieval of genomic regions)
    # If your input file is not sorted by chromosome and position, which is required by tabix1. You can check if your file is sorted by using the sort command with the -c option:
#############################################################

zcat mm10.refGene.gtf.gz | sort -c -k1,1 -k4,4n -k5,5n -t$'\t'
    
    # If the output says sort: disorder: …, then you need to sort your file by chromosome and position before compressing and indexing it with tabix2. You can do this by using the following commands:

zcat mm10.refGene.gtf.gz | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip > mm10.refGene.gtf.gz.sorted.gtf.gz

tabix -p gff mm10.refGene.gtf.gz.sorted.gtf.gz

#########################################

zcat mm10.ncbiRefSeq.gtf.gz | sort -c -k1,1 -k4,4n -k5,5n -t$'\t'
    
    # If the output says sort: disorder: …, then you need to sort your file by chromosome and position before compressing and indexing it with tabix2. You can do this by using the following commands:

zcat mm10.ncbiRefSeq.gtf.gz | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip > mm10.ncbiRefSeq.sorted.gtf.gz
tabix -p gff mm10.ncbiRefSeq.sorted.gtf.gz

#########################################

zcat mm10.ensGene.gtf.gz | sort -c -k1,1 -k4,4n -k5,5n -t$'\t'
    
    # If the output says sort: disorder: …, then you need to sort your file by chromosome and position before compressing and indexing it with tabix2. You can do this by using the following commands:

zcat mm10.ensGene.gtf.gz | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip > mm10.ensGene.sorted.gtf.gz
tabix -p gff mm10.ensGene.sorted.gtf.gz

#########################################
zcat mm10.knownGene.gtf.gz | sort -c -k1,1 -k4,4n -k5,5n -t$'\t'
    
    # If the output says sort: disorder: …, then you need to sort your file by chromosome and position before compressing and indexing it with tabix2. You can do this by using the following commands:

zcat mm10.knownGene.gtf.gz | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip > mm10.knownGene.sorted.gtf.gz
tabix -p gff mm10.knownGene.sorted.gtf.gz

# WGS samples; for the region of Pcdh_abg clusters; chr18:36726000-38419200

# Basic plot by using only bam files
############################################################
    modbamtools plot -r chr18:36726000-38419200 \
    --gtf gencode.vM22.annotation.sorted.gtf.gz \
    --out . \
    --prefix WGS_sample_01_WT_Pcdh_abg \
    --samples WGS_sample_01_WT \
    --track-titles Genes\
    WT.merged.sup.mod.sort.bam 
############################################################

############################################################
    modbamtools plot -r chr18:36726000-38419200 \
    --gtf gencode.vM22.annotation.sorted.gtf.gz \
    --out . \
    --prefix WGS_sample_02_Hom_Pcdh_abg \
    --samples WGS_sample_02_Hom \
    --track-titles Genes\
    Hom.merged.sup.mod.sort.bam 
############################################################


# Adding Bigwig tracks; Bigwig tracks can also be added to the plot with --bigwig:
############################################################
    modbamtools plot -r chr18:36726000-38419200 \
    --gtf gencode.vM22.annotation.sorted.gtf.gz \
    --out . \
    --bigwig H3K4me3WT1_rpkm.bigWig \
    --bigwig H3K4me3Miss1_rpkm.bigWig \
    --prefix WGS_sample_01_WT_Pcdh_abg_H3K4me3 \
    --samples WGS_sample_01_WT \
    --track-titles Genes,H3K4me3_WT1,H3K4me3_Miss1\
    WT.merged.sup.mod.sort.bam 
############################################################
    modbamtools plot -r chr18:36726000-38419200 \
    --gtf gencode.vM22.annotation.sorted.gtf.gz \
    --out . \
    --bigwig H3K4me3WT1_rpkm.bigWig \
    --bigwig H3K4me3Miss1_rpkm.bigWig \
    --prefix WGS_sample_02_Hom_Pcdh_abg_H3K4me3 \
    --samples WGS_sample_01_Hom \
    --track-titles Genes,H3K4me3_WT1,H3K4me3_Miss1\
    Hom.merged.sup.mod.sort.bam 
############################################################


#  We can add more samples to the plot
############################################################
    modbamtools plot -r chr18:36726000-38419200 \
    --gtf gencode.vM22.annotation.sorted.gtf.gz \
    --out . \
    --bigwig H3K4me3WT1_rpkm.bigWig \
    --bigwig H3K4me3Miss1_rpkm.bigWig \
    --prefix WT-WGS_Hom-WGS_Pcdh_abg_H3K4me3 \
    --samples WGS_sample_01_WT,WGS_sample_02_Hom \
    --track-titles Genes,H3K4me3_WT1,H3K4me3_Miss1\
    WT.merged.sup.mod.sort.bam\
    Hom.merged.sup.mod.sort.bam\
############################################################        


