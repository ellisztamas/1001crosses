#!/usr/bin/env bash

# Calculate SNP heritabilities using GEMMA.
#
# Tom Ellis, 16th September 2025

# SLURM
#SBATCH --job-name=01_calculate_h2
#SBATCH --output=05_results/06_SNP_heritability/slurm/%x-%a.out
#SBATCH --error=05_results/06_SNP_heritability/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-11

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh



# ===Input files ==== #

i=$SLURM_ARRAY_TASK_ID

# Array of text files with Blup phenotypes.
blup_array=(03_processing/06_process_phenotypes/output/*blups*tsv)
pheno_file=${blup_array[$i]}

# Array of VCF files for the two F8 replicates and the parents
parents_vcf=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
progeny_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz


# === Output files ===

# Output directory
outdir=$scratchdir/06_SNP_heritability
mkdir -p $outdir

# Basename for the file containing the PVE estimate.
pve_file=$(basename -s .tsv $pheno_file)_pve.txt

# === Script ===

if [[ "$pheno_file" == *"parent"* ]]; then
    echo "Phenotype filename contains 'parent'"
    input_vcf=$parents_vcf
else
    echo "Phenotype filename does not contain 'parent'"
    input_vcf=$progeny_vcf
fi

echo "VCF file: ${input_vcf}"
echo "Phenotype file: ${pheno_file}"
echo "Output directory: ${outdir}"
echo ""

echo "Calculating SNP heritability."
02_library/run_SNP_heritability.sh \
  --vcf $input_vcf \
  --phenotypes $pheno_file \
  --outdir $outdir

cp $outdir/$pve_file 05_results/06_SNP_heritability/output/