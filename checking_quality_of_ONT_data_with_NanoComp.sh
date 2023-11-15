# Quality control of ONT data, including read length, number of total reads, number of sequenced base etc...

NanoComp --summary \
a01_sequencing_summary.txt \
a02_sequencing_summary.txt \
a03_sequencing_summary.txt \

--names a01 a02 a03 --outdir compare-runs --threads 128
