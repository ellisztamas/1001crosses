#!/usr/bin/env bash

# Merge VCF files from individual samples into a single VCF file.
# Rename samples based on IBDpainting results where necesary.

# SLURM
#SBATCH --mem=2GB
#SBATCH --job-name=08_merge_VCF
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --qos=short
#SBATCH --time=2:00:00



# Set working directory and load conda environment
source setup.sh

# === Input ===

# Directory with VCF for individual samples.
indir=$scratchdir/04_resequencing/05_snp_calls

# Validation results.
ibdpainting=03_processing/04_resequencing/output/ibdpainting_results_resequenced.tsv


# === Output ===

# Output directory
outdir=$scratchdir/04_resequencing/08_merge_VCF
mkdir -p $outdir

# Sample list for correcting names.
# This lists old names, and those to replace, separated by a space
sample_list=$outdir/sample_list.txt

# VCF file with merged samples, and incorrect names
vcf_merged=$outdir/merged_samples.vcf.gz
# VCF file with merged samples, and correct names
outfile=$outdir/resequenced.vcf.gz


# === Main ===

# Sample list for correcting names.
grep -e "correct\|swap\|missing haplotype\|ambiguous Swede" $ibdpainting | \
awk -F'\t' '{print $4, " ", $5}' \
> $sample_list

# Merge VCF files that start with a number (i.e. aren't blank)
bcftools merge \
    --output-type u \
    --output $vcf_merged \
    $indir/[0-9]*.vcf.gz

# Swap sample names in the merged VCF file
bcftools reheader \
    -s $sample_list \
    -o $outfile \
    $vcf_merged
tabix $outfile

cp ${outfile}* 03_processing/04_resequencing/output/