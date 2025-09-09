#!/usr/bin/env bash

# Run GWAS on simulated datasets
# 
# The previous scipt creates many independent simulated datasets for parents,
# both F8 cohorts, and permuted datasets. This script runs GEMMA on each of
# them.
#
# Tom Ellis, 21st July 2025

# SLURM
#SBATCH --job-name=05_run_gemma
#SBATCH --output=05_results/22_simulated_datasets/slurm/%x-%a.out
#SBATCH --error=05_results/22_simulated_datasets/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=60:00
#SBATCH --array=0-3200
# Set working directory and load the conda environment
source setup.sh


# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# Array of phenotype files.
phenotype_file_array=($scratchdir/22_simulated_datasets/04_simulate_phenotypes/*phen)
phenotype_file=${phenotype_file_array[$i]}

# Array of paths to VCF files.
# For permuted datasets there is only a directory path - the file names will be
# added later
genotype_array=(03_processing/05_imputation/output/parental_lines.vcf.gz
    03_processing/05_imputation/output/F8_cohort1_phased_imputed.vcf.gz
    03_processing/05_imputation/output/F8_cohort2_phased_imputed.vcf.gz
    ${scratchdir}/22_simulated_datasets/02_permute_parents
    )

# === Output ===

outdir=$scratchdir/22_simulated_datasets/${SLURM_JOB_NAME}
mkdir -p $outdir

# Additional arguments passed to GEMMA
# Run all three kinds of statistical test using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"

# === Main === 

# Which VCF file to use depends on whether the phenotype file starts with 
# `parents`, `cohort1`, `cohort2` or `permutd`
filename=$(basename ${phenotype_file})
case "$filename" in
    parents*)
        echo "Processing using parental genotypes"
        vcf_file=${genotype_array[0]}
        ;;
    cohort1*)
        echo "Processing using genotypes for cohort1"
        vcf_file=${genotype_array[1]}
        ;;
    cohort2*)
        echo "Processing using genotypes for cohort2"
        vcf_file=${genotype_array[2]}
        ;;
    permutd*)
        dataset=$(echo $filename | cut -d"_" -f1,2)
        echo "Processing using genotypes for permuted dataset ${dataset}"
        vcf_file=${genotype_array[3]}/${dataset}.vcf
        ;;
    *)
        echo "Filename does not match of the expected strings"
        ;;
esac

# Run GWAS
02_library/run_GEMMA_with_without_K.sh \
  --vcf $vcf_file \
  --phenotypes $phenotype_file \
  --outdir $outdir \
  --gemma_args "${gemma_args}"
