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
conda activate gemma

# Paths to input files
vcf=03_processing/pieters_sample_sheet/output/F8_snp_matrix.vcf.gz
pheno_file=05_results/03_gwas_ft10/input/ft10_rep1.tsv
covariates=05_results/03_gwas_ft10/input/covariates_rep1.tsv
outdir=05_results/03_gwas_ft10/output
# Additional arguments passed to GEMMA
gemma_args="-lmm 4 -maf 0.05"

02_library/run_GEMMA.sh \
  --vcf $vcf \
  --phenotypes $pheno_file \
  --outdir $outdir \
  --covariates $covariates \
  --gemma_args "${gemma_args}"

date