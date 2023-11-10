# To execute the contents of the .bashrc file in the current shell session (necessary to start conda)
    source ~/.bashrc

# To create bam.bai files (Generic indexer for TAB-delimited genome position files) install Samptools, bcftools and HTSlib 
# If you are working on Shark CentOS Slurm cluster, use shark terminal window to load these modules

module load library/htslib/1.16/gcc-8.5.0
module load genomics/ngs/bcftools/1.11/gcc-8.3.1
module load genomics/ngs/samtools/1.16.1/gcc-8.5.0

# To check installed modules
module li

# This will take the first million reads and give you counts of each read length. To look at more of the bam, get rid of the "head" command
# This command outputs the length distribution of both mapped and unmapped reads
samtools view E4.merged.sup.mod.sort.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n | uniq -c

# If you want to save the output of all reads and their lengths to a text file, you can modify your command to process all reads and redirect the output to a .txt file. 
samtools view E4.merged.sup.mod.sort.bam | \
cut -f 10 | \
perl -ne 'chomp;print length($_) . "\n"' | \
sort -n | \
uniq -c > reads_lengths.txt

# The samtools view command without any flags will output both mapped and unmapped reads from the BAM file. 
# If you want to get the length distribution of only the mapped reads, you need to use the -F 4 flag with samtools view, which excludes unmapped reads.

samtools view -F 4 E4.merged.sup.mod.sort.bam | \
cut -f 10 | \
perl -ne 'chomp;print length($_) . "\n"' | \
sort -n | \
uniq -c > mapped_reads_lengths.txt
