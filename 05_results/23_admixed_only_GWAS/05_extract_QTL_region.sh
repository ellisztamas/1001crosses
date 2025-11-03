#!/usr/bin/env bash

# Create a text file giving genotypes around the QTL for seed size on chr 5.
#
# Tom Ellis, 15th October 2025

# Load conda environment
source setup.sh 

# === Input files === #

# Directory with flowering-time BLUPs
phenotype_file=05_results/23_admixed_only_GWAS/output/01_subset_lines/seed_size_blups_NxS1.tsv
phenotype_name=${phenotype_name}

# VCF files for parents and F8s
input_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz


# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=05_results/23_admixed_only_GWAS/output/05_extract_QTL_region
mkdir -p $outdir

genotype_array=$outdir/Chr5_QTL_region.txt


# === Main === #

sample_names=$(cut -f2 $phenotype_file | paste -sd, -)

echo "snp,"${sample_names} > $genotype_array

bcftools view \
    --samples $sample_names \
    --regions Chr5:16207115-16248614 \
    -i "MAC >=10" \
    $input_vcf |
bcftools query \
        --format "%CHROM\_%POS[,%GT]\n" |
    sed 's/0\/0/0/g; s/0|0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/0|1/1/g; s/1|0/1/g; s/1\/1/2/g; s/1|1/2/g; s/\.\/./NA/g' \
    >> $genotype_array