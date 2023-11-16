
# Nanopore Data Processing and Visualization for bam files

## Overview
This repository focuses on processing and visualizing Oxford nanopore long-read sequencing data, particularly for detecting DNA methylation. It includes tools such as Methylartist and Modbamtools for parsing, manipulating, and visualizing methylation patterns and DNA/RNA base modification data stored in BAM format.

Methylartist: tools for visualizing modified bases from nanopore sequence data
Modbamtools is a set of tools to manipulate and visualize DNA/RNA base modification data that are stored in bam format

## Contents
1. **CpG_density_calculation_for_genes.r**: An R script for calculating CpG density in genes.
2. **Intro_to_ONT_sequencing.pptx**: A PowerPoint presentation introducing Oxford Nanopore Technologies (ONT) sequencing.
3. **bam_files_filtering_more_than_150bp_length_reads.sh**: A shell script for filtering BAM files to retain reads longer than 150bp.
4. **bam_files_to_fastq_files.sh**: Script to convert BAM files to FASTQ format.
5. **bam_files_to_fastq_files_with_alignment_to_reference.sh**: Script to convert BAM files to FASTQ format with alignment to a reference genome.
6. **checking_quality_of_ONT_data_with_NanoComp.sh**: A script for checking the quality of ONT data using NanoComp.
7. **filtering_mapping_quality.sh**: Script for filtering based on mapping quality.
8. **merging_BAM_files_with_samtools.sh**: Script for merging BAM files using samtools.
9. **methylartist_all_commands_for_bam.sh**: A comprehensive shell script for using Methylartist with BAM files.
10. **methylartist_for_bam_files.sh**: A script for using Methylartist specifically with BAM files.
11. **modbam2bed_running_sample.sh**: Script for running Modbamtools sample.
12. **modbamtools_for_bam_files.sh**: A script for using Modbamtools with BAM files.
13. **number_of_CpG_loci_in_GRCm38_genome.R**: An R script to count the number of CpG loci in the GRCm38 mouse genome.
14. **number_of_reads_with_length.sh**: Script to count the number of reads with a specified length.

## Installation and Usage
Details on installation and usage are not explicitly provided in the repository. Users are likely expected to have a working knowledge of shell scripting and R programming.

## References
- Methylartist: [GitHub](https://github.com/adamewing/methylartist), [Publication](https://doi.org/10.1093/bioinformatics/btac292)
- Modbamtools: [GitHub](https://rrazaghi.github.io/modbamtools/), [Preprint](https://www.biorxiv.org/content/10.1101/2022.07.07.499188v1)

---

Please note: This README is based on the file names and types in the repository and may not fully represent the intended use or capabilities of the scripts and tools provided. For detailed information, please refer to the individual files.
