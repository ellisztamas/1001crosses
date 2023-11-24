# Code and links to download the raw data files

# load packages
module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate 1001crosses

# A. thaliana TAIR10 reference genome
tair10=https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
wget -P $01_data/01_reference_genome $tair10
# index the genome now
bwa index $genome
