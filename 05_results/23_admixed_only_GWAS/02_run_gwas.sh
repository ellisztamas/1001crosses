#!/usr/bin/env bash

# Run GEMMA on subsets of F9 lines that are descended from Northern and Southern
# parents.
#
# Tom Ellis, 13th October 2025

# SLURM
#SBATCH --job-name=02_run_gwas
#SBATCH --output=05_results/23_admixed_only_GWAS/slurm/%x-%a.out
#SBATCH --error=05_results/23_admixed_only_GWAS/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=0-5

# Load conda environment
source setup.sh 

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory with flowering-time BLUPs
pheno_file_array=(05_results/23_admixed_only_GWAS/output/01_subset_lines/*tsv)
phenotype_file=${pheno_file_array[$i]}
phenotype_name=${phenotype_name}

# VCF files for parents and F8s
input_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz


# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=05_results/23_admixed_only_GWAS/output/02_run_gwas
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
  --input $outdir/with_K/${phenotype_name}.assoc.txt \
  --outDir $outdir/with_K
02_library/plot_gwas.py \
  --input $outdir/no_K/${phenotype_name}.assoc.txt \
  --outDir $outdir/no_K

# Manhattan plots as an interactive HTML, which are slow
02_library/interactive_manhattan_plot.R \
    --input $outdir/with_K/${phenotype_name}.assoc.txt \
    --threshold 2
02_library/interactive_manhattan_plot.R \
    --input $outdir/no_K/${phenotype_name}.assoc.txt \
    --threshold 2
