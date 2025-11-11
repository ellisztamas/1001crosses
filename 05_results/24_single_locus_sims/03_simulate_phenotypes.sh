#!/usr/bin/env bash

# Apply the R script to simulate phenotypes to genotype matrices.
# 
# The R script 02_simulate_phenotypes.R simulates phenotypes based on genotype-
# matrix files for several levels of heritability. This applies that script
# to genotype matrix files for different levels of heritability.
# 
# Tom Ellis, 10th November 2025

# SLURM
#SBATCH --job-name=03_simulate_phenotypes
#SBATCH --output=05_results/24_single_locus_sims/slurm/%x.out
#SBATCH --error=05_results/24_single_locus_sims/slurm/%x.err
#SBATCH --qos=rapid
#SBATCH --time=10:00

# Set working directory and load the conda environment
source setup.sh


# === Input ===

# Array of genotype-matrix files for the parents, for different MAF bins.
array_of_parental_files=(05_results/24_single_locus_sims/output/01_samples_genotypes/parents_genotype_matrix_*.tsv)

# The genotype matrix files have a total of N loci (rows), for k replicate datasets.
# Each simulation should only use n loci, so n*k must equal N
# In this case n=1, so this is easy
nloci=1


# === Output ===

outdir=05_results/24_single_locus_sims/output/03_simulate_phenotypes
mkdir -p $outdir



# === Main ===

# There are four files each for each cohort, corresponding to allele frequency bins.
for i in 0 1 2 3
do
    # File paths for each cohort for a single MAF bin.
    parents=${array_of_parental_files[$i]}
    cohort1=${parents/parents/cohort1}
    cohort2=${parents/parents/cohort2}
    # Run the Rscript.
    Rscript 05_results/24_single_locus_sims/02_simulate_phenotypes.R \
        --parents $parents \
        --cohort1 $cohort1 \
        --cohort2 $cohort2 \
        --nloci $nloci \
        --output $outdir
done