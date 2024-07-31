#!/usr/bin/env bash

# Concatenate summary results from the simulations.
#
# Inputs:
#     Files giving details and summary results for each GEMMA run
# Output:
#     A single file with this information in one place

# Tom Ellis, 12th July 2024

# SLURM
#SBATCH --job-name=concat_sims
#SBATCH --output=05_results/12_correlated_background/slurm/%x.out
#SBATCH --error=05_results/12_correlated_background/slurm/%x.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=20:00

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

# === Input === #

# indir=05_results/12_correlated_background/output/tmp
indir=/scratch-cbe/users/thomas.ellis/crosses/05_results/12_correlated_background
file_array=($indir/*/*/*/*/sim_summary.txt)

# === Output === #

outfile=05_results/12_correlated_background/output/simulation_summary.csv

# === Script ===

echo "chr,pos,cohort,liability,includes_K,p_target_SNP,beta,gif,n_false_pos" > $outfile

for f in ${file_array[@]}
do
    cat $f >> $outfile
done

