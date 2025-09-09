#!/usr/bin/env bash

# Run the R script to merge result files
#
# Tom Ellis, 23rd June 2025

# SLURM
#SBATCH --job-name=merge_gemma_results
#SBATCH --output=05_results/19_ft12_resampled/slurm/%x.out
#SBATCH --error=05_results/19_ft12_resampled/slurm/%x.err
#SBATCH --mem=20GB
#SBATCH --qos=short
#SBATCH --time=4:00:00

# Set working directory and load the conda environment
source setup.sh

mkdir -p 05_results/19_ft12_resampled/output
Rscript 05_results/19_ft12_resampled/03_join_results.R