#!/usr/bin/env bash

# PCA on the two replicates of F8s plus the parents.
#
# Use PLINK to prune and filter the VCF files for the parents and both replicates
# of F8s, then run a PCA on each.
# 
# Tom Ellis, adapting code from  https://speciationgenomics.github.io/pca/, 
# 15th December 2023

# SLURM
#SBATCH --job-name=pca
#SBATCH --output=05_results/01_pca/slurm/%x-%a.out
#SBATCH --error=05_results/01_pca/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-1

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

i=$SLURM_ARRAY_TASK_ID

# ===Input files ==== #

# VCF files for parents and F8s
parents=03_processing/05_imputation/output/parental_lines.vcf.gz
progeny=03_processing/05_imputation/output/F8_phased_imputed.vcf.gz
# Make a Bash array of the three VCF files
vcf_files=($parents $progeny)
# File to use in this job
infile=${vcf_files[$i]}



# === Output files ===

# Output directory
outdir=05_results/01_pca/output
mkdir -p $outdir
# Suffix for output files
file_suffix=$outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )



# === Script ===

# perform linkage pruning - i.e. identify prune sites
# 
# window size of 50kb, window step size of 10bp, pruning SNPs with r2 > 0.2
# Select loci with minor-allele frequency > 5%, with no more than 10 missing data
# Select individuals with no more than 50% missing data
plink \
  --vcf $infile \
  --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --maf 0.05 --geno 0.1 --mind 0.5 \
  --indep-pairwise 50 10 0.4 \
  --out $file_suffix

# Run the PCA
plink \
    --vcf $infile \
    --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.05 --geno 0.1 --mind 0.5 \
    --pca 'header' \
    --extract $file_suffix.prune.in \
    --out $file_suffix

# Tidy up the output files.
rm $file_suffix.{prune.in,prune.out,nosex}