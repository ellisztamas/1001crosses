#!/usr/bin/env bash

# Script to run 01_simulate_phenotypes.py as a SLURM job

# Tom Ellis, 15th July 2024

# SLURM
#SBATCH --job-name=sim_phenotypes
#SBATCH --output=05_results/12_correlated_background/slurm/%x.out
#SBATCH --error=05_results/12_correlated_background/slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

# === Input === #

# indir=05_results/11_sim_background/output/tmp
indir=/scratch-cbe/users/thomas.ellis/crosses/05_results/12_correlated_background

# === Output === #

pheno_file_list=$indir/pheno_file_list.txt

# === Script === #

# Run the script
python 05_results/12_correlated_background/01_simulate_phenotypes.py

# Save a list of output files for later, because this takes ages.
find $indir/**/parents/**/phenotype_file.tsv > $pheno_file_list