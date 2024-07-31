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

# Tom Ellis, 10th July 2024

# SLURM
#SBATCH --job-name=sim_gemma
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-299

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
phenotype_dir=(05_results/09_gwas_single_loci/output/tmp/*)
snp_dir=${phenotype_dir[$i]}

# Additional arguments passed to GEMMA
# Run using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"

# === Output files ===

# Output directory for results
parents_dir=$snp_dir/parents
rep1_dir=$snp_dir/rep1
rep2_dir=$snp_dir/rep2

# === Script === 

echo "GWAS on the parents."
02_library/run_GEMMA_with_without_K.sh \
  --vcf $parental_snp_matrix \
  --phenotypes $parents_dir/phenotype_file.tsv \
  --outdir $parents_dir \
  --gemma_args "${gemma_args}"
  # rm -r $parents_dir/tmp
# Create some plots
# 02_library/plot_gwas.py \
#   --input $parents_dir/with_K/phenotype_file.assoc.txt \
#   --outDir $parents_dir/with_K

# 02_library/plot_gwas.py \
#   --input $parents_dir/no_K/phenotype_file.assoc.txt \
#   --outDir $parents_dir/no_K


echo "GWAS on replicate 1 of the F8s"
02_library/run_GEMMA_with_without_K.sh \
  --vcf $F8_rep1_snp_matrix \
  --phenotypes $rep1_dir/phenotype_file.tsv \
  --outdir $rep1_dir \
  --gemma_args "${gemma_args}"
# rm -r $rep1_dir/tmp

# 02_library/plot_gwas.py \
#     --input $rep1_dir/with_K/phenotype_file.assoc.txt \
#     --outDir $rep1_dir/with_K
# 02_library/plot_gwas.py \
#   --input $rep1_dir/no_K/phenotype_file.assoc.txt \
#   --outDir $rep1_dir/no_K


echo "GWAS on replicate 2 of the F8s"
02_library/run_GEMMA_with_without_K.sh \
  --vcf $F8_rep2_snp_matrix \
  --phenotypes $rep2_dir/phenotype_file.tsv \
  --outdir $rep2_dir \
  --gemma_args "${gemma_args}"
# rm -r $rep2_dir/tmp
# Create some plots
# 02_library/plot_gwas.py \
#   --input $rep1_dir/with_K/phenotype_file.assoc.txt \
#   --outDir $rep1_dir/with_K
# 02_library/plot_gwas.py \
#   --input $rep1_dir/no_K/phenotype_file.assoc.txt \
#   --outDir $rep1_dir/no_K