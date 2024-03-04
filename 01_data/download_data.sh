# Code and links to download the raw data files

# load conda environment
source setup.sh

# GATK needs to be loaded locally
ml build-env/f2021
ml gatk/4.2.6.1-gcccore-10.2.0-java-13

# A. thaliana TAIR10 reference genome
tair10=https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
wget -P 01_data/01_reference_genome $tair10
# index the genome now
bwa index $genome
gatk CreateSequenceDictionary -R 01_data/01_reference_genome/TAIR10_chr_all.fas