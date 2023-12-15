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
source 03_processing/pieters_sample_sheet/00_setup.sh

# Input data
# Directory containing SNP calls for each chromosome separately.
indir=$workdir/05_snp_calls
# VCF files for each chromosome separately
VCF_chr1=$indir/F8_snp_matrix_chr1.vcf.gz
VCF_chr2=$indir/F8_snp_matrix_chr2.vcf.gz
VCF_chr3=$indir/F8_snp_matrix_chr3.vcf.gz
VCF_chr4=$indir/F8_snp_matrix_chr4.vcf.gz
VCF_chr5=$indir/F8_snp_matrix_chr5.vcf.gz

# Target files
# File name for the concatenated SNP matrix
vcf_with_full_paths=$indir/F8_snp_matrix_full_paths.vcf.gz
vcf_with_sample_names=$indir/F8_snp_matrix.vcf.gz
# Directory to stage out the resulting SNP matrix
outdir=03_processing/pieters_sample_sheet/output

# Concatenate VCF files
bcftools concat \
    -o $vcf_with_full_paths \
    $VCF_chr1 $VCF_chr2 $VCF_chr3 $VCF_chr4 $VCF_chr5
if [ $? -eq 0 ] ; then echo "Concatenating VCF files completed successfully"; fi

# Rename the samples to remove absolute paths and .bam suffixes
bcftools query -l $vcf_with_full_paths > $indir/header_names.txt
xargs -rd '\n' -a $indir/header_names.txt basename -a --suffix=.bam > $indir/new_header_names.txt
bcftools reheader \
    -s $indir/new_header_names.txt \
    -o $vcf_with_sample_names \
    $vcf_with_full_paths
if [ $? -eq 0 ] ; then echo "Renamed samples successfully"; fi
# Index
tabix $vcf_with_sample_names

# Stage out
cp $vcf_with_sample_names $outdir
cp $vcf_with_sample_names.tbi $outdir