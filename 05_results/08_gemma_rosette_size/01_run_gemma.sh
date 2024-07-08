#!/usr/bin/env bash

# Run a vanilla univariate GWA using GEMMA on rosette size with and without a
# kinship matrix
#
# Data are from measurements at 49 days from F9 and parental plants grown in a 
# randomised experiment at 12°C
#
# Inputs:
#     VCF file for the parents
#     VCF file for the F8s
#     Text files giving phenotype data on rosette size showing:
#         A column of 1s to tell GEMMA to fit an intercept.
#         A column of genotype IDs
#         A column of phenotypes
# Output:
#     Output files from GEMMA for each replicate

# Tom Ellis, 25th June 2024

# SLURM
#SBATCH --job-name=rosette_size_gemma
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --array=0-2

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

# === Paths to input files ===

i=$SLURM_ARRAY_TASK_ID

# VCF files for parents and F8s
F8_rep1_snp_matrix=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.vcf.gz
F8_rep2_snp_matrix=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.vcf.gz
parental_snp_matrix=03_processing/04_pieters_VCF/output/parental_snp_matrix.vcf.gz
# Array of paths to VCF files.
# It is really important that these are in the same order as the phenotype files!
vcf_array=($F8_rep1_snp_matrix $F8_rep2_snp_matrix $parental_snp_matrix)
input_vcf=${vcf_array[$i]}

# Array of phenotype files for the two sets of replicate crosses plus the parents. 
gemma_input_files=(03_processing/06_process_phenotypes/output/rosette_size_*.tsv)
in_phenotype=${gemma_input_files[$i]}

# Additional arguments passed to GEMMA
# Run all three kinds of statistical test using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"

# === Output files ===

# Output directory for GEMMA results and temporary files.
outdir=05_results/08_gemma_rosette_size/output
mkdir -p $outdir

# === Script === 

echo "\nVCF file: ${input_vcf}"
echo "Phenotype file: ${in_phenotype}"
echo "Output directory: ${outdir}\n\n"

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

date