#!/usr/bin/env bash

# Create intermediate files for simulation and GWAS.
# 
# Script to convert VCF files for observed (non-permuted) VCF files for the 
# parents and F8s to Plink format.
# (For permuted datasets, this is done on creation in 02_permute_parents.sh.)
#
# Tom Ellis, 21st July 2025

# SLURM
#SBATCH --job-name=03_create_plink
#SBATCH --output=05_results/22_simulated_datasets/slurm/%x.out
#SBATCH --error=05_results/22_simulated_datasets/slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=1:00:00

# Set working directory and load the conda environment
source setup.sh


# === Input ===

# VCF files
# Parents
parents_vcf='03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz'
# Two cohorts of F8s
cohort1_vcf='03_processing/09_impute_haplotypes/output/F8_cohort1_phased_imputed.vcf.gz'
cohort2_vcf='03_processing/09_impute_haplotypes/output/F8_cohort2_phased_imputed.vcf.gz'


# === Output ===

outdir=$scratchdir/22_simulated_datasets/${SLURM_JOB_NAME}
mkdir -p $outdir

# Prefix for output files without any file extension
parents_out=$outdir/$(basename -s .vcf.gz ${parents_vcf})
cohort1_out=$outdir/$(basename -s .vcf.gz ${cohort1_vcf})
cohort2_out=$outdir/$(basename -s .vcf.gz ${cohort2_vcf})

tmp_snplist=$outdir/tmp_snplist.txt

# === Main === 

# Create Plink files
echo "Converting parents to Plink format."
plink2 \
  --vcf $parents_vcf \
  --maf 0.01 \
  --make-bed \
  --set-all-var-ids '@_#' \
  --out $parents_out

# Extract variant IDs (2nd column of the .bim) that survived the filter
awk '{print $2}' "${parents_out}.bim" > "$tmp_snplist"
echo "  -> $(wc -l < "$tmp_snplist") loci kept"

# Convert the F8s, using only the SNPs that vary in the parents.
echo "Converting cohort 1 to Plink format."
plink2 \
  --vcf $cohort1_vcf \
  --extract $tmp_snplist \
  --make-bed \
  --set-all-var-ids '@_#' \
  --out $cohort1_out

echo "Converting cohort 2 to Plink format."
plink2 \
  --vcf $cohort2_vcf \
  --extract $tmp_snplist \
  --make-bed \
  --set-all-var-ids '@_#' \
  --out $cohort2_out

