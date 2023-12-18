#!/usr/bin/env bash

# PCA on the two replicates of F8s plus the parents.
#
# Use PLINK to prune and filter the VCF files for the parents and both replicates
# of F8s, then run a PCA on each.
# 
# Input: Three VCF files
# Output: Eigenvector, eigenvalue and log files for each VCF.
#
# Tom Ellis, adapting code from  https://speciationgenomics.github.io/pca/, 
# 15th December 2023

# SLURM
#SBATCH --job-name=pca
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-2

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate 1001crosses

i=$SLURM_ARRAY_TASK_ID

# Path to a VCF files with dubious samples removed
rep1=03_processing/pieters_VCF/F8_snp_matrix_purged_rep1.vcf.gz
rep2=03_processing/pieters_VCF/F8_snp_matrix_purged_rep2.vcf.gz
# VCF file for the parents
parents=01_data/03_parental_genotypes/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz
# Make a Bash array of the three VCF files
vcf_files=($rep1 $rep2 $parents)
# File to use in this job
infile=${vcf_files[$i]}

# Output directory
outdir=05_results/01_pca/output
mkdir -p $outdir
# Suffix for output files
file_suffix=$outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )

# perform linkage pruning - i.e. identify prune sites
# 
# window size of 50kb, window step size of 10bp, pruning SNPs with r2 > 0.4
# Select loci with minor-allele frequency > 5%, with no more than 10$ missing data
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