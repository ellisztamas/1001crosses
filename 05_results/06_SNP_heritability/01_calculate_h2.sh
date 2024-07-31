#!/usr/bin/env bash

# Calculate SNP heritabilities using GEMMA.
#
# Tom Ellis, 28th June 2024

# SLURM
#SBATCH --job-name=snp_h2
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=2

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

i=$SLURM_ARRAY_TASK_ID

# ===Input files ==== #

# Array of VCF files for the two F8 replicates and the parents
vcf_dir=03_processing/04_pieters_VCF/output
vcf_files=($vcf_dir/F8_snp_matrix_purged_rep1.vcf.gz $vcf_dir/F8_snp_matrix_purged_rep2.vcf.gz 03_processing/04_pieters_VCF/output/parental_snp_matrix.vcf.gz)
# File to use in this job
input_vcf=${vcf_files[$i]}

# Array of files containing genetic values for rosette size
blup_dir=03_processing/06_process_phenotypes/output
rs_files=($blup_dir/rosette_size_blups_F9_rep1.tsv $blup_dir/rosette_size_blups_F9_rep2.tsv $blup_dir/rosette_size_blups_parents.tsv)
input_rs=${rs_files[$i]}

# Array of files containing genetic values for seed size
ss_files=($blup_dir/seed_size_blups_F9_rep1.tsv $blup_dir/seed_size_blups_F9_rep2.tsv $blup_dir/seed_size_blups_parents.tsv)
input_ss=${ss_files[$i]}

# === Output files ===

# Output directory
outdir=05_results/06_SNP_heritability/output
mkdir -p $outdir

# === Script ===

echo "\nVCF file: ${input_vcf}"
echo "Rosette size phenotype file: ${input_rs}"
echo "Output directory: ${outdir}\n\n"

echo "Calculating SNP heritability for rosette size"
02_library/run_SNP_heritability.sh \
  --vcf $input_vcf \
  --phenotypes $input_rs \
  --outdir $outdir

echo "SNP heritability for seed size"
02_library/run_SNP_heritability.sh \
  --vcf $input_vcf \
  --phenotypes $input_ss \
  --outdir $outdir

date