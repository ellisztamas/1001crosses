#!/usr/bin/env bash

# Run GEMMA on simulated datasets with and without the correction for the 
# relatedness matrix.
#
# Inputs:
#     VCF file for the parents
#     VCF file for the F8s
#     Text files with simulated phenotypes from 04_simulate_phenotypes.py
# Output:
#     Output files from GEMMA for each replicate

# Tom Ellis, 12th July 2024

# SLURM
#SBATCH --job-name=sim_corr_background
#SBATCH --output=05_results/12_correlated_background/slurm/%x-%a.out
#SBATCH --error=05_results/12_correlated_background/slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=0-1199

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

# Array of phenotype files for the two sets of replicate crosses plus the parents. 
file_array=( $(cat /scratch-cbe/users/thomas.ellis/crosses/05_results/12_correlated_background/pheno_file_list.txt) )
# A single file name for the parents and F8s
parent_file=${file_array[$i]}
rep1_file=${parent_file/parents/rep1}
rep2_file=${parent_file/parents/rep2}


# Additional arguments passed to GEMMA
# Run using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"

# === Output files === #

# Output directories
outdir_parents=$(dirname $parent_file)
outdir_rep1=$(dirname $rep1_file)
outdir_rep2=$(dirname $rep2_file)

# === Script === 

echo "GWAS on the parents."
02_library/run_GEMMA_with_without_K.sh \
  --vcf $parental_snp_matrix \
  --phenotypes $parent_file \
  --outdir $outdir_parents \
  --gemma_args "${gemma_args}"


echo "GWAS on replicate 1 of the F8s"
02_library/run_GEMMA_with_without_K.sh \
  --vcf $F8_rep1_snp_matrix \
  --phenotypes $rep1_file \
  --outdir $(dirname $rep1_file) \
  --gemma_args "${gemma_args}"


echo "GWAS on replicate 2 of the F8s"
02_library/run_GEMMA_with_without_K.sh \
  --vcf $F8_rep2_snp_matrix \
  --phenotypes $rep2_file \
  --outdir $(dirname $rep2_file) \
  --gemma_args "${gemma_args}"