#!/usr/bin/env bash

# SNP calling, based on a SNP matrix created by Fernando Rabanal.


# 

# Input: aligned BAM files
# Output: five zipped VCF files.

# Tom Ellis, 5th March 2024, adapting code from this tutorial by Kathryn Hodgins:
# https://khodgins.github.io/Bioinformatics_Introduction/Topic_5/

# SLURM
#SBATCH --job-name=call_SNPs_gatk
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=40GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00

date

# Setup
# Terminate the script automatically if any part of it fails (returns a non-zero exit status)
set -e
# Set the working directory and activate the conda environment
source setup.sh
# Load the modules for GATK manually
ml build-env/f2021
ml gatk/4.2.6.1-gcccore-10.2.0-java-13

# Inputs
# Directory containing aligned BAM files
indir=$workdir/05_base_recalibration
# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas
# SNP matrix file
parental_SNP_matrix=03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix_mac20.vcf.gz

# Outputs
# Output directory
outdir=$workdir/05_snp_calls_GATK
mkdir -p $outdir
# file with location of bamfiles (one per line)
bam_list=${outdir}/samplelist.txt

# Script
# create file with bamfile paths
ls -d ${indir}/[1-9]*x*.bam > $bam_list

# Run HaplotypeCaller on each sample
for infile in `cat ${bam_list} | head -n 1`
do
outfile=$(basename ${infile/.bam/.g.vcf})
gatk --java-options "-Xmx15g" HaplotypeCaller \
   -R $genome \
   -I $infile \
   --native-pair-hmm-threads 3 \
   -ERC GVCF \
   -O $outdir/$outfile
done
if [ $? -eq 0 ] ; echo "HaplotypeCaller completer successfully" ; fi