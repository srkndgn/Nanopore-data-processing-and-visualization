#!/bin/bash
#SBATCH --job-name=methylartist
#SBATCH --time=96:0:0 # hh:mm:ss

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=160000 # 160G
#SBATCH --mail-user=s.dogan@lumc.nl
#SBATCH --mail-type=ALL
#
#SBATCH --comment=Devepi

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!

__conda_setup="$('/exports/humgen/Serkan/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/path/to/conda_directory/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/path/to/conda_directory/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/path/to/conda_directory/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup

# <<< conda initialize <<<

conda activate methylartist


set -e # exit run_tests.sh if anything fails


# locus
# Pcdh_a_b_g clusters

time methylartist locus -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -i chr18:36771904-37909338 -r GRCm38.p6.genome.fa -g gencode.vM22.annotation.sorted.gtf.gz -n CG  -p 3,5,2,3,3 --width 30 --height 16 --svg
echo -e "\done...\n"
time methylartist locus -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -i chr18:36771904-37909338 -r GRCm38.p6.genome.fa -g gencode.vM22.annotation.sorted.gtf.gz -n CG  -p 3,5,2,3,3 --width 30 --height 16
echo -e "\done...\n"


# region
# Pcdh_a_b_g clusters

# all with genes names
time methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --labelgenes --skip_align_plot --panelratio 7,0,3,7 --svg
echo -e "\done...\n"

time methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --labelgenes --skip_align_plot --panelratio 7,0,3,7
echo -e "\done...\n"

# all without genes names
time methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --skip_align_plot --panelratio 7,0,3,7 --svg
echo -e "\done...\n"
time methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a8522 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --skip_align_plot --panelratio 7,0,3,7
echo -e "\done...\n"

# just genes with names
time methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --labelgenes --skip_align_plot --panelratio 10,0,1,0 --svg
echo -e "\done...\n"
time methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --labelgenes --skip_align_plot --panelratio 10,0,1,0
echo -e "\done...\n"

# just genes without names
time methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --skip_align_plot --panelratio 10,0,1,0 --svg
echo -e "\done...\n"
time methylartist region -i chr18:36771904-37909338 -b sample_01.sorted.bam:#9a32a8,sample_02.sorted.bam:#32a852 -g gencode.vM22.annotation.sorted.gtf.gz -p 32 -n CG -r GRCm38.p6.genome.fa --skip_align_plot --panelratio 10,0,1,0
echo -e "\done...\n"


