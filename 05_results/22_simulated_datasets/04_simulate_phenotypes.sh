#!/usr/bin/env bash

# Simulate phenotype files for GWAS
# 
# Simulate data files with normally-distributed phenotypes for the parents, both
# cohorts of F8s, and replicate permuted datasets using GCTA.
# Phenotypes have
#     heritabilities of 0.2, 0.4, 0.6 or 0.8
#     1, 5, 50 or 500 causal loci
#     200 replicates

# For each replicate the initial simulation is done for the parents, then that
# list of causal SNPs and effect sizes is used for the other datasets so that the
# results can be directly compared.
#
# Tom Ellis, 22nd July 2025

# SLURM
#SBATCH --job-name=04_simulate_phenotypes
#SBATCH --output=05_results/22_simulated_datasets/slurm/%x-%a.out
#SBATCH --error=05_results/22_simulated_datasets/slurm/%x-%a.err
#SBATCH --mem=5 GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=1-50

# Set working directory and load the conda environment
source setup.sh


# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# Genotype files.
# Variables are given in Plink format without a file extension
# Parents
parents_in=${scratchdir}/22_simulated_datasets/03_create_plink/parental_lines
# Two cohorts of F8s
cohort1_in=${scratchdir}/22_simulated_datasets/03_create_plink/F8_cohort1_phased_imputed
cohort2_in=${scratchdir}/22_simulated_datasets/03_create_plink/F8_cohort2_phased_imputed
# A single replicate of permuted parental genotypes.
permutd_in=${scratchdir}/22_simulated_datasets/02_permute_parents/permutd_$(printf "%03d\n" $i)
# An array you can loop over
genotype_array=($parents_in $cohort1_in $cohort2_in $permutd_in)


# === Output ===

outdir=$scratchdir/22_simulated_datasets/${SLURM_JOB_NAME}
mkdir -p $outdir

# Text file contaiing all SNPs in the parental .bim file
all_snps=$outdir/all_snps.txt

# Addition intermediate files are defined in the loops below, because they
# depend on loop parameters:
    # SNP markers to use
    # Effect sizes at each SNP
    # Phenotype file

phenotype_array=("parents" "cohort1" "cohort2" "permutd")

# === Main ===

# A list of all SNPs in the parental genotype file.
awk '{print $2}' "${parents_in}.bim" > $all_snps

# Loop over values for number of loci and heritability, and simulate phenotypes.
for nsnps in 1 5 50 500; do
    echo "Simulating phenotypes with $nsnps SNPs."
    # Randomly select 50 SNPs
    selected_snps=${outdir}/snps_${i}_${nsnps}.txt
    shuf $all_snps | head -$nsnps > $selected_snps
    
    # Assign random effect sizes (normal distribution, mean=0, sd=1)
    effect_sizes=${outdir}/effect_sizes_${i}_${nsnps}.txt
    awk '{printf "%s\t%.4f\n", $1, (rand()*2-1)}' $selected_snps > $effect_sizes

    # Simulate phenotypes for different levels of heritability
    for heritability in 0.2 0.4 0.6 0.8; do
        echo "\tSimulating phenotypes with heritability ${heritability}."
        for base in {0..3}; do
            pheno_out=$outdir/${phenotype_array[$base]}_$(printf "%03d\n" $i)_${nsnps}_$heritability
            gcta64 \
                --bfile "${genotype_array[$base]}" \
                --simu-qt \
                --simu-causal-loci $effect_sizes \
                --simu-hsq $heritability \
                --out $pheno_out 
        done
    done
done
