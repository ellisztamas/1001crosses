#!/usr/bin/env bash

# Run GEMMA on resampled F9 flowering-time datasets
#
# Tom Ellis, 23rd June 2025

# SLURM
#SBATCH --job-name=19_ft12_resampled
#SBATCH --output=05_results/19_ft12_resampled/slurm/%x-%a.out
#SBATCH --error=05_results/19_ft12_resampled/slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-199

# Load conda environment
source setup.sh 

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory with data
# There is one working directory for each resampled dataset.
dir_array=($scratchdir/19_ft12_resampled/ft12_resample_*)
indir=${dir_array[$i]}

# File with resampled F9 phenotypes
phenotype_file=${indir}/phenotypes.tsv
# VCF files all F8s
input_vcf=03_processing/05_imputation/output/F8_phased_imputed.vcf.gz


# === Output files === #

# Save results to the same place as the phenotype file.
outdir=$indir

# Additional arguments passed to GEMMA
# Run all three kinds of statistical test using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"



# === Main === #


echo "VCF file: ${input_vcf}"
echo "Phenotype file: ${phenotype_file}"
echo "Output directory: ${outdir}"

# GWAS with the kinship matrix
02_library/run_GEMMA_with_without_K.sh \
  --vcf $input_vcf \
  --phenotypes $phenotype_file \
  --outdir $outdir \
  --gemma_args "${gemma_args}"
