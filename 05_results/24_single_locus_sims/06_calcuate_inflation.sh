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
#SBATCH --job-name=06_calculate_inflation
#SBATCH --output=05_results/24_single_locus_sims/slurm/%x.out
#SBATCH --error=05_results/24_single_locus_sims/slurm/%x.err
#SBATCH --mem=48G
#SBATCH --qos=medium
#SBATCH --time=24:00:00


# Set working directory and load the conda environment
source setup.sh

# === Inputs === 

# Array of results files
noK_indir=$scratchdir/24_single_locus_sims/04_run_gemma/no_K
withK_indir=$scratchdir/24_single_locus_sims/04_run_gemma/with_K

# Path to the script to calculate lambda
lambda=05_results/24_single_locus_sims/05_calculate_inflation.R


# === Outputs ===

outdir=scratchdir/24_single_locus_sims/06_calculate_inflation
mkdir -p $outdir

# Final file to write to
noK_outfile=05_results/24_single_locus_sims/output/inflation_noK.tsv
withK_outfile=05_results/24_single_locus_sims/output/inflation_withK.tsv


# === Main === 


# File headers giving cohort, replicate ID, minor-allele freq range, heritability,
# number of significant windows, and lambda statistic, all tab separated.
echo -e "cohort\trep\tmaf_range\th2\tlambda\tnwindows" > $noK_outfile
echo -e "cohort\trep\tmaf_range\th2\tlambda\tnwindows" > $withK_outfile

# Run the script on each results file
# Results without a K matrix
for infile in ${noK_indir}/*assoc.txt; do
    echo "Processing file $(basename ${infile})"
    Rscript $lambda --in $infile --out $noK_outfile
done
# Results with a K matrix
for infile in ${withK_indir}/*assoc.txt; do
    echo "Processing file $(basename ${infile})"
    Rscript $lambda --in $infile --out $withK_outfile
done
