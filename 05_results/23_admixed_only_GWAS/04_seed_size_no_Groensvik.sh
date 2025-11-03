#!/usr/bin/env bash

# Run GEMMA on subsets of F9 lines that are descended from Northern and Southern
# parents, but excluding lines with a parent from the Grönsvik peninsula.
#
# Tom Ellis, 13th October 2025

# SLURM
#SBATCH --job-name=04_seed_size_no_Groensvik
#SBATCH --output=05_results/23_admixed_only_GWAS/slurm/%x-%a.out
#SBATCH --error=05_results/23_admixed_only_GWAS/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-1

# Load conda environment
source setup.sh 

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory with phenotype files
pheno_file_array=(05_results/23_admixed_only_GWAS/output/01_subset_lines/seed_size*tsv)
phenotype_file=${pheno_file_array[$i]}
phenotype_name=$(basename -s '.tsv' $phenotype_file)

# VCF files for parents and F8s
input_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz


# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=05_results/23_admixed_only_GWAS/output/04_seed_size_no_Groensvik
mkdir -p $outdir

# Smaller phenotype file excluding lines with a parent from Groensvik
subsetted_phenotype_file=$outdir/${phenotype_name}_subset.tsv

# Additional arguments passed to GEMMA
# Run all three kinds of statistical test using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"


# === Main === #


echo "VCF file: ${input_vcf}"
echo "Phenotype file: ${phenotype_file}"
echo "Output directory: ${outdir}"

# Remove the lines from Grönsvik from the phenotype file.
groensvik_lines="1435x6240_rep1\|6020x6030_rep1\|6020x6030_rep2\|6220x9339_rep1\|6220x9339_rep2\|9388x9057_rep1\|9388x9057_rep2"
grep -v $groensvik_lines $phenotype_file > $subsetted_phenotype_file

# GWAS with the kinship matrix
02_library/run_GEMMA_with_without_K.sh \
  --vcf $input_vcf \
  --phenotypes $subsetted_phenotype_file \
  --outdir $outdir \
  --gemma_args "${gemma_args}"

# Create some plots
# Manhattan and QQ plots as .png files, which are quick
02_library/plot_gwas.py \
  --input $outdir/with_K/${phenotype_name}_subset.assoc.txt \
  --outDir $outdir/with_K
02_library/plot_gwas.py \
  --input $outdir/no_K/${phenotype_name}_subset.assoc.txt \
  --outDir $outdir/no_K

# Manhattan plots as an interactive HTML, which are slow
02_library/interactive_manhattan_plot.R \
    --input $outdir/with_K/${phenotype_name}_subset.assoc.txt \
    --threshold 2
02_library/interactive_manhattan_plot.R \
    --input $outdir/no_K/${phenotype_name}_subset.assoc.txt \
    --threshold 2
