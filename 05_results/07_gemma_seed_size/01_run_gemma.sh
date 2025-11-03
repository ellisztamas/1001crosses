#!/usr/bin/env bash

# Run a vanilla univariate GWA using GEMMA on seed-size BLUPs.
#
# Tom Ellis, 11th September 2025

# SLURM
#SBATCH --job-name=01_run_gemma
#SBATCH --output=05_results/07_gemma_seed_size/slurm/%x-%a.out
#SBATCH --error=05_results/07_gemma_seed_size/slurm/%x-%a.err
#SBATCH --mem=20GB
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=0-3

# Load conda environment
source setup.sh 

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory with flowering-time BLUPs
pheno_file_array=(03_processing/06_process_phenotypes/output/seed_size*tsv)
phenotype_file=${pheno_file_array[$i]}

# VCF files for parents and F8s
progeny_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz
parents_vcf=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
# Array of paths to VCF files.
# It is really important that these are in the same order as the phenotype files!
vcf_array=($progeny_vcf $parents_vcf $progeny_vcf $progeny_vcf)
input_vcf=${vcf_array[$i]}



# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=05_results/07_gemma_seed_size/output
mkdir -p $outdir
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

# Create some plots
# Manhattan and QQ plots as .png files, which are quick
02_library/plot_gwas.py \
  --input $outdir/with_K/$(basename -s'.tsv' ${phenotype_file}).assoc.txt \
  --outDir $outdir/with_K
02_library/plot_gwas.py \
  --input $outdir/no_K/$(basename -s'.tsv' ${phenotype_file}).assoc.txt \
  --outDir $outdir/no_K

# Manhattan plots as an interactive HTML, which are slow
02_library/interactive_manhattan_plot.R \
    --input $outdir/with_K/$(basename -s'.tsv' ${phenotype_file}).assoc.txt \
    --threshold 2
02_library/interactive_manhattan_plot.R \
    --input $outdir/no_K/$(basename -s'.tsv' ${phenotype_file}).assoc.txt \
    --threshold 5
