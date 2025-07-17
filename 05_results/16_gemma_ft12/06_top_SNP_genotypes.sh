#!/usr/bin/env bash

# Pull out the genotypes for flowering time SNPs.
# 
# This takes a file of candidate SNPs and pulls out genotypes at each.
# Note that the output has SNPs as rows and samples as columns, and you probably
# want the transpose of that. Heterozygotes are converted to NA.
#
# Tom Ellis, 1st July 2025

# SLURM
#SBATCH --job-name=06_top_SNP_genotypes
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


snp_list_csv=05_results/16_gemma_ft12/output/candidate_peak_positions.csv


# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=$scratchdir/16_gemma_ft12
mkdir -p $outdir

# The list of SNPs, but with tabs as the delimiter, for BCFtools
snp_list_tsv=$outdir/top_SNPs/candidate_peak_positions.tsv

# Output files from VCF tools.
# This has a row for each marker and a column for each SNP
# Genotypes are in VCF format (0|0, 0|1, ./, etc)
vcftools_GT_format=$outdir/$(basename -s .vcf.gz $input_vcf)_top_SNPS
# The same data, but with genotypes converted to integers (or NA)
integer_genotypes=${vcftools_GT_format}.tsv

# === Main ===

# Convert SNP list to be tab-delimited
awk -F',' '{print "Chr"$1, $2}' OFS='\t' $snp_list_csv > $snp_list_tsv

# Pull genotypes at focal loci in VCF format
vcftools \
    --gzvcf $input_vcf \
    --positions $snp_list_tsv \
    --extract-FORMAT-info GT \
    --out $vcftools_GT_format

# Convert VCF genotypes to integers.
cat ${vcftools_GT_format}.GT.FORMAT | \
    sed \
        -e 's_0|0_0_g' \
        -e 's_0|1_1_g' \
        -e 's_1|0_1_g' \
        -e 's_1|1_2_g' \
        -e 's_./._NA_g' \
    > $integer_genotypes

cp $integer_genotypes 05_results/16_gemma_ft12/output/top_SNPs