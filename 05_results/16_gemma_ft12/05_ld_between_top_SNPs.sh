#!/usr/bin/env bash

# Calculate LD between SNPs associated with flowering time.
#
# This takes a hand-curated list of SNPs that look like they are at the top of a
# peak for flowering time, and get LD between all of them.
#
# Tom Ellis, 30th June 2025

# SLURM
#SBATCH --job-name=05_ld_between_top_SNPs
#SBATCH --output=05_results/16_gemma_ft12/slurm/%x-%a.out
#SBATCH --error=05_results/16_gemma_ft12/slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-3

# Set working directory and load the conda environment
source setup.sh


# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# VCF files for parents and F8s
progeny_vcf=03_processing/05_imputation/output/F8_phased_imputed.vcf.gz
parents_vcf=03_processing/05_imputation/output/parental_lines.vcf.gz
cohort1_vcf=03_processing/05_imputation/output/F8_cohort1_phased_imputed.vcf.gz
cohort2_vcf=03_processing/05_imputation/output/F8_cohort2_phased_imputed.vcf.gz
# Array of paths to VCF files.
# It is really important that these are in the same order as the phenotype files!
vcf_array=($progeny_vcf $parents_vcf $cohort1_vcf $cohort2_vcf)
input_vcf=${vcf_array[$i]}


snp_list=05_results/16_gemma_ft12/output/top_SNPs/candidate_peak_positions.csv


# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=$scratchdir/16_gemma_ft12
mkdir -p $outdir

plink_prefix=$outdir/$(basename -s .vcf.gz $input_vcf)

# Final LD file. This is tab-delimited
ld_file=05_results/16_gemma_ft12/output/top_SNPs/$(basename -s .vcf.gz $input_vcf).ld


# === Main ===

# Filter SNPs and set SNP labels to match the CSV file (chr,ps)
plink2 \
    --vcf $input_vcf \
    --set-all-var-ids '@,#' \
    --extract $snp_list \
    --make-bed \
    --out $plink_prefix
# Create a table of LD between all pairs
plink \
    --bfile $plink_prefix \
    --r2 inter-chr \
    --ld-window-r2 0 \
    --out $plink_prefix

# Substitute whitespace for tabs
awk '{$1=$1}1' OFS='\t' ${plink_prefix}.ld > $ld_file
