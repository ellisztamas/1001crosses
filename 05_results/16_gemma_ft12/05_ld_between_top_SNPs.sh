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
#SBATCH --mem=4GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-3

# Set working directory and load the conda environment
source setup.sh


# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory with flowering-time BLUPs
pheno_file_array=(03_processing/06_process_phenotypes/output/flowering_time*tsv)
phenotype_file=${pheno_file_array[$i]}
pheno_basename=$(basename -s .tsv $phenotype_file)

# VCF files for parents and F8s
progeny_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz
parents_vcf=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
# Array of paths to VCF files.
# It is really important that these are in the same order as the phenotype files!
vcf_array=($progeny_vcf $parents_vcf $progeny_vcf $progeny_vcf)
input_vcf=${vcf_array[$i]}

# A hand-made list of markers that look like the tips of peaks
snp_list=05_results/16_gemma_ft12/output/top_SNPs/candidate_peak_positions.csv


# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=$scratchdir/16_gemma_ft12/05_ld_between_top_SNPs
mkdir -p $outdir

# Prefic for intermediate files from plink
plink_prefix=$outdir/${pheno_basename}

# Final LD file. This is tab-delimited
mkdir -p 05_results/16_gemma_ft12/output/05_ld_between_top_SNPs
ld_file=05_results/16_gemma_ft12/output/05_ld_between_top_SNPs/${pheno_basename}.ld


# === Main ===

echo "Filtering individuals and markers to match those used in the GWAS."
plink2 \
    --vcf $input_vcf \
    --set-all-var-ids '@,#' \
    --extract $snp_list \
    --make-bed \
    --max-alleles 2 \
    --keep $phenotype_file \
    --out $plink_prefix

echo "Creating a table of LD between all pairs."
plink \
    --bfile $plink_prefix \
    --r2 inter-chr \
    --ld-window-r2 0 \
    --out $plink_prefix

# Substitute whitespace for tabs
awk '{$1=$1}1' OFS='\t' ${plink_prefix}.ld > $ld_file
