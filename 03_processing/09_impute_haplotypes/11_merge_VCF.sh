#!/usr/bin/env bash

# Merge spliced VCF files for individual F8s into a single VCF file.
#
# This also writes a list of line names that is used to ensure phenotype files
# match.
#
# Tom Ellis, 11th August 2025

# SLURM
#SBATCH --job-name=11_merge_VCF
#SBATCH --output=03_processing/09_impute_haplotypes/slurm/%x.out
#SBATCH --error=03_processing/09_impute_haplotypes/slurm/%x.err
#SBATCH --mem=5GB
#SBATCH --qos=medium
#SBATCH --time=24:00:00

# Set working directory
source setup.sh
# set -e

# === Input ===

# Directory containing individual VCF files for each F8 individual that have
# been imputed from parental genotypes and recombination breakpoints.
indir=$scratchdir/09_impute_haplotypes/10_splice_VCF


# === Output ===

outdir=03_processing/09_impute_haplotypes/output
mkdir -p $outdir

# Temporary text file to store input VCF file names.
vcf_file_list=$outdir/vcf_file_list.txt

# VCF containing all merged individuals
merged_vcf=$outdir/F8_imputed.vcf.gz
# List of line names in the output file.
sample_names=${merged_vcf/.vcf.gz/line_names.txt}

# === Main ===

echo "Merging VCF files."
bcftools merge \
    --output-type z \
    --output $merged_vcf \
    ${indir}/*_sorted_haplotypes.vcf.gz

echo "Writing a list of line names."
bcftools query -l $merged_vcf > $sample_names