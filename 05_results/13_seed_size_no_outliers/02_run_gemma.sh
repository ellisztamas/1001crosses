#!/usr/bin/env bash

# Run univariate GWA using GEMMA on seed size with and without a
# kinship matrix of relatedness. Data are excluding lines with unusually large
# seeds. See 01_identify_outliers.R.
# 
#
# Inputs:
#     VCF file for the parents
#     VCF files for the two replicates of F8s
#     Text files giving phenotype data on seed size, excluding outliers, showing:
#         A column of 1s to tell GEMMA to fit an intercept.
#         A column of genotype IDs
#         A column of phenotypes
# Output:
#     Output files from GEMMA for each replicate

# Tom Ellis, 17th July 2024

# SLURM
#SBATCH --job-name=seed_size_no_outliers
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=0-2

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

# === Paths to input files ===

i=$SLURM_ARRAY_TASK_ID

# VCF files for parents and F8s
parental_snp_matrix=03_processing/04_pieters_VCF/output/parental_snp_matrix.vcf.gz
F8_rep1_snp_matrix=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.vcf.gz
F8_rep2_snp_matrix=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.vcf.gz
# Array of paths to VCF files.
# It is really important that these are in the same order as the phenotype files!
vcf_array=($parental_snp_matrix $F8_rep1_snp_matrix $F8_rep2_snp_matrix)
input_vcf=${vcf_array[$i]}

# Array of phenotype files for the two sets of replicate crosses plus the parents. 
gemma_input_files=(05_results/13_seed_size_no_outliers/output/seed_size_*.tsv)
in_phenotype=${gemma_input_files[$i]}

# Additional arguments passed to GEMMA
# Run all three kinds of statistical test using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"

# === Output files ===

# Output directory for results
outdir=05_results/13_seed_size_no_outliers/output
mkdir -p $outdir

# === Script === 

echo "\nVCF file: ${input_vcf}"
echo "Phenotype file: ${in_phenotype}"
echo "Output directory: ${outdir}"

# GWAS with the kinship matrix
02_library/run_GEMMA_with_without_K.sh \
  --vcf $input_vcf \
  --phenotypes $in_phenotype \
  --outdir $outdir \
  --gemma_args "${gemma_args}"

# Create some plots
02_library/plot_gwas.py \
  --input $outdir/with_K/$(basename -s'.tsv' ${in_phenotype}).assoc.txt \
  --outDir $outdir/with_K

02_library/plot_gwas.py \
  --input $outdir/no_K/$(basename -s'.tsv' ${in_phenotype}).assoc.txt \
  --outDir $outdir/no_K