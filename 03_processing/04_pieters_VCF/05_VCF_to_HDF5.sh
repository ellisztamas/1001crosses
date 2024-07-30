#!/usr/bin/env bash

# Convert VCF files for the parents and replicate cohorts of F8s to HDF5 using
# scikit allel.
#
# Inputs:
#     VCF file for the parents
#     VCF files for both cohorts of F8s
# Output:
#     HDF5 files for each VCF

# Tom Ellis, 5th July 2024

# SLURM
#SBATCH --job-name=vcf_to_hdf5
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --array=0-1

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
vcf_array=($F8_rep1_snp_matrix $F8_rep2_snp_matrix $parental_snp_matrix)
input_vcf=${vcf_array[$i]}

# === Output files ===

outfile=${input_vcf/.vcf.gz/.hdf5}

# === Script === 

python 02_library/vcf_to_HDF5.py \
    --input $input_vcf --output $outfile