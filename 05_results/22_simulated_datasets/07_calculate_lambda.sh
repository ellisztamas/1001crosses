#!/usr/bin/env bash

# Calculate lambda statistic for each GWAS result file
# 
# Run an Rscript to calculate lambda on each result file.
# Each script outputs a row giving cohort ID (parents, cohorts 1 or 2, or 
# permutd), replicate ID, number of loci, heritability and lambda.
# This row is appended to a text file.
#
# Tom Ellis, 21st July 2025

# SLURM
#SBATCH --job-name=07_calculate_lambda
#SBATCH --output=05_results/22_simulated_datasets/slurm/%x.out
#SBATCH --error=05_results/22_simulated_datasets/slurm/%x.err
#SBATCH --mem=48G
#SBATCH --qos=medium
#SBATCH --time=24:00:00


# Set working directory and load the conda environment
source setup.sh

# === Inputs === 

# Array of results files for the parents only
input_dir=$scratchdir/22_simulated_datasets/05_run_gemma/no_K

# Path to the script to calculate lambda
lambda=05_results/22_simulated_datasets/06_calculate_lambda.R

# === Outputs ===

outdir=scratchdir/22_simulated_datasets/07_calculate_lambda
mkdir -p $outdir


# Final file to write to
outfile=05_results/22_simulated_datasets/output/lambda_values_noK.tsv


# === Main === 


# File headers giving the locus names and a column for r2
echo -e "cohort\trep\tnloci\th2\tlambda" > $outfile

# Run the script on each results file
for infile in ${input_dir}/*assoc.txt; do
    echo "Processing file $(basename ${infile})"
    Rscript $lambda --in $infile --out $outfile
done
