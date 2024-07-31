#!/usr/bin/env bash

# Script to run 01_allele_correlations.py as a SLURM job

# Tom Ellis, 8th July 2024

# SLURM
#SBATCH --job-name=af_correlations
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

python 05_results/09_gwas_single_loci/01_allele_correlations.py