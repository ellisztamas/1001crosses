#!/usr/bin/env bash

# Use PLINK to calculate relatedness matrices between the parents and F8s.
#
# Tom Ellis, 22nd May 2025

# SLURM
#SBATCH --job-name=01_calculate_matrices
#SBATCH --output=05_results/18_relatedness_matrix/slurm/%x-%a.out
#SBATCH --error=05_results/18_relatedness_matrix/slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=rapid
#SBATCH --time=10:00
#SBATCH --array=0-1

date 

source setup.sh

# === Input files ===

i=$SLURM_ARRAY_TASK_ID

# VCF files for parents and F8s
parents=03_processing/05_imputation/output/parental_lines.vcf.gz
progeny=03_processing/05_imputation/output/F8_phased_imputed.vcf.gz
# Make a Bash array of the three VCF files
vcf_files=($parents $progeny)
# File to use in this job
infile=${vcf_files[$i]}


# === Output === 

outdir=05_results/18_relatedness_matrix/output
mkdir -p $outdir

# Filename prefix for the output files
out_prefix=$outdir/$(basename -s .vcf.gz $infile)


# === Main === 

# Create a triangular matrix of relatedness
# Filtering for MAF > 0.05, missingness per SNP < 10%, and missingness per individual < 50%
# Use the --make-rel square0 option to create a square matrix
# Use the --cov option to skip standardization by allele frequency
plink \
    --vcf $infile \
    --maf 0.05 --geno 0.1 --mind 0.5 \
    --make-rel square0 cov \
    --out $out_prefix