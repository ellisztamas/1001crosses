#!/usr/bin/env bash

# Calculate allele frequencies and Hardy-Weinberg statistics.

# Input:
#     VCF files for the parents, and replicates one and two of the offspring
# Output:
#     .frq: Minor-allele frequency report

# Tom Ellis, 18th December 2023

# SLURM
#SBATCH --job-name=allele_freqs
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-2

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Path to a VCF files with dubious samples removed
rep1=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.vcf.gz
rep2=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.vcf.gz
# VCF file for the parents
parents=03_processing/04_pieters_VCF/output/parental_snp_matrix.vcf.gz
# Make a Bash array of the three VCF files
vcf_files=($rep1 $rep2 $parents)
# File to use in this job
infile=${vcf_files[$i]}

# === Output files === #

# Output directory
outdir=05_results/05_allele_freqs/output
mkdir -p $outdir
# Suffix for Plink output files
file_suffix=$outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )



# === Script === #

plink2 \
  --vcf $infile \
  --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --freq \
  --out $file_suffix

vcftools --gzvcf $infile --freq2 --out $file_suffix

date

