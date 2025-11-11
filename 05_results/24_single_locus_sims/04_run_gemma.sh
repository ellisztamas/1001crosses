#!/usr/bin/env bash

# 
# This assumes that files have the same samples!
#
# Tom Ellis, 10th November 2025

# SLURM
#SBATCH --job-name=04_run_gemma
#SBATCH --output=05_results/24_single_locus_sims/slurm/%x-%a.out
#SBATCH --error=05_results/24_single_locus_sims/slurm/%x-%a.err
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=2-1199

# Set working directory and load the conda environment
source setup.sh


# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# Array of phenotype files with heritability = 0.2
dataset_array=(05_results/24_single_locus_sims/output/03_simulate_phenotypes/*0.2.tsv)
# Matching input files with different levels of heritability
pheno_02=${dataset_array[$i]}
pheno_04=${pheno_02/0.2./0.4.}
pheno_06=${pheno_02/0.2./0.6.}
pheno_08=${pheno_02/0.2./0.8.}

# Set an input VCF file to use depending on whether the phenotypes are for 
# parental or F8 genotypes.
if [[ "$pheno_02" == *"parent"* ]]; then
  vcf_file=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
elif [[ "$pheno_02" == *"cohort"* ]]; then
  vcf_file=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz
else
  echo "Error: pheno_02 ('$pheno_02') does not contain 'parents' or 'cohort'." >&2
  exit 1 # Terminate the script with an error status
fi
echo "vcf_file is set to: $vcf_file"

# Additional arguments passed to GEMMA
# Run all three kinds of statistical test using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"

# === Output ===

outdir=$scratchdir/24_single_locus_sims/04_run_gemma
mkdir -p $outdir


# === Main ===

heritability_array=($pheno_02 $pheno_04 $pheno_06 $pheno_08)

for phenotype_file in ${heritability_array[@]}
do
    echo "Running GEMMA with file ${phenotype_file}."
    02_library/run_GEMMA_with_without_K.sh \
        --vcf $vcf_file \
        --phenotypes $phenotype_file \
        --outdir $outdir \
        --gemma_args "${gemma_args}"
done