#!/usr/bin/env bash

# Combine per chromosome VCF files into one.
# Rename the samples in the VCF files to remove absolute paths and suffixes

# Input: VCF files for each chromosome separately for all the samples
# Output: A single VCF file

# Tom Ellis, adapting code by Pieter Clauw, 27th November 2023

# SLURM
#SBATCH --job-name=merge_VCF
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=01:00:00

# Set working directory
source setup.sh

# === Input data === #

# Directory containing SNP calls for each chromosome separately.
indir=$scratchdir/07_snp_calls
# VCF files for each chromosome separately
VCF_chr1=$indir/F8_snp_matrix_Chr1.vcf.gz
VCF_chr2=$indir/F8_snp_matrix_Chr2.vcf.gz
VCF_chr3=$indir/F8_snp_matrix_Chr3.vcf.gz
VCF_chr4=$indir/F8_snp_matrix_Chr4.vcf.gz
VCF_chr5=$indir/F8_snp_matrix_Chr5.vcf.gz

# === Define outputs === #

# Directory for the output files
outdir=$scratchdir/08_merge_VCF
mkdir -p $outdir
# File name for the concatenated SNP matrix
vcf_with_full_paths=$outdir/F8_snp_matrix_full_paths.vcf.gz
vcf_with_sample_names=$outdir/F8_snp_matrix.vcf.gz
# Directory to stage out the resulting SNP matrix
projdir=03_processing/02_original_sample_sheet/output

# === Script === #

# Concatenate VCF files
rm $vcf_with_full_paths
bcftools concat \
    -o $vcf_with_full_paths \
    $VCF_chr1 $VCF_chr2 $VCF_chr3 $VCF_chr4 $VCF_chr5
if [ $? -eq 0 ] ; then echo "Concatenating VCF files completed successfully"; fi

# Rename the samples to remove absolute paths and .bam suffixes
# Create text files giving the old and new names
bcftools query -l $vcf_with_full_paths > $outdir/header_names.txt # Old names
python 03_processing/02_original_sample_sheet/07_rename_samples.py # New names
# Swap the names
rm $vcf_with_sample_names
bcftools reheader \
    -s $outdir/new_header_names.txt \
    -o $vcf_with_sample_names \
    $vcf_with_full_paths
if [ $? -eq 0 ] ; then echo "Renamed samples successfully"; fi
# Index the VCF file
tabix $vcf_with_sample_names

# Stage out
cp $vcf_with_sample_names $projdir
cp $vcf_with_sample_names.tbi $projdir