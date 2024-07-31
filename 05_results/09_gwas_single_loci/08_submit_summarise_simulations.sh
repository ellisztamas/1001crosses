#!/usr/bin/env bash

# Script to run 07_summarise_sims.py as a SLURM job

# Tom Ellis, 11th July 2024

# SLURM
#SBATCH --job-name=summarise_sims
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=8:00:00

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

python 05_results/09_gwas_single_loci/07_summarise_sims.py