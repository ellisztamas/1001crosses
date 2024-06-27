#!/usr/bin/env bash

# Run a vanilla univariate GWA using GEMMA on flowering time phenotypes at 10°C.
#
# Tom Ellis, 13th December 2023

# SLURM
#SBATCH --job-name=ft10_gwas
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=8:00:00

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate 1001crosses

# Paths to input files
F8_snp_matrix=03_processing/02_original_sample_sheet/output/F8_snp_matrix.vcf.gz
parental_snp_matrix=03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix_mac20.vcf.gz
gemma_input_files=05_results/04_gwas_ft10/gemma_input_files
# Output directory for GEMMA results and temporary files.
outdir=05_results/04_gwas_ft10/output
mkdir -p $outdir
# Additional arguments passed to GEMMA
# Run all three kinds of statistical test using a minor-allele-frequency of 0.05
gemma_args="-lmm 4 -maf 0.05"

echo "Running GEMMA on replicate 1 of the F8s..."
02_library/run_GEMMA.sh \
  --vcf $F8_snp_matrix \
  --phenotypes $gemma_input_files/ft10_rep1.tsv \
  --outdir $outdir \
  --covariates $gemma_input_files/ft10_rep1_covariates.tsv \
  --gemma_args "${gemma_args}"
if [ $? -eq 0 ] ; then echo "done.\n\n"; fi

date

echo "Running GEMMA on replicate 2 of the F8s..."
02_library/run_GEMMA.sh \
  --vcf $F8_snp_matrix \
  --phenotypes $gemma_input_files/ft10_rep2.tsv \
  --outdir $outdir \
  --covariates $gemma_input_files/ft10_rep2_covariates.tsv \
  --gemma_args "${gemma_args}"
if [ $? -eq 0 ] ; then echo "done.\n\n"; fi

date

echo "Running GEMMA on the parents..."
02_library/run_GEMMA.sh \
  --vcf $parental_snp_matrix \
  --phenotypes $gemma_input_files/ft10_parents.tsv \
  --outdir $outdir \
  --gemma_args "${gemma_args}"
if [ $? -eq 0 ] ; then echo "done.\n\n"; fi

date