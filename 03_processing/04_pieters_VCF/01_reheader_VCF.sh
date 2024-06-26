#!/usr/bin/env bash

# Reheader the VCF for the F8s with sample names and split by replicate.

# Pieter's VCF file uses absolute paths to raw BAM files in the header
# This script pulls them out in order, replaces them with sample names, and 
# creates a new VCF file these shorter headers. Based on these names it then 
# creates separate VCF files for F8 replicates 1 and 2.

# Inputs:
#     Pieter's VCF file for the F8 plants
#     Pieter's calls on which samples could be validated (indirectly: this is done
#         in an R script, which this script calls)
# Outputs:
#     Zipped VCF files for the two replicates, excluding samples that could not
#         be validated, and with sample names as headers.

# Tom Ellis, 15th December 2023

# SLURM
#SBATCH --job-name=process_pieters_VCF
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=20GB
#SBATCH --qos=rapid
#SBATCH --time=30:00

date 

source ./setup.sh

# == Inputs ==

# Path to the VCF file for the F8s created by Pieter
input_vcf=../crosses_continued/004.F8/001.genotyping/003.results/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC_crossesF8.vcf.gz

# === Output files ===

# Filename for the old headers
old_header=03_processing/04_pieters_VCF/output/vcf_header_to_change.txt
# Filename to save the new headers. This is created inside the R script, so don't mess with it.
new_header=03_processing/04_pieters_VCF/output/new_vcf_header.txt
# File with samples to be kept
samples_to_keep=03_processing/04_pieters_VCF/output/samples_to_keep.txt
# Path to save the resulting VCF file.
vcf_with_all_samples=03_processing/04_pieters_VCF/output/F8_snp_matrix.vcf.gz
# Path to save a VCF file with dubious samples removed
purged_VCF=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged.vcf.gz

# === Script ===

# Pull out the existing header
bcftools query -l $input_vcf > $old_header
# Replace absolute paths with sample names, and save as $new_header
# Conda does not play nicely with r-tidyverse, so you may need to run this separately.
Rscript 03_processing/04_pieters_VCF/02_reorder_header.R

# # Swap the headers and keep only the samples confirmed by SNPmatch
echo "Swapping headers and removing dubious genotypes."
bcftools reheader --samples $new_header $input_vcf > $vcf_with_all_samples
bcftools view --samples-file $samples_to_keep $vcf_with_all_samples > $purged_VCF

# Bash witchcraft to generate a comma-separated list of samples for each replicate:
# 1. extract sample names,
# 2. select only samples with 'rep1' in the name
# 3. Convert from having one sample per row to a comma-separated list.
# 4. Remove the trailing comma.
echo "Splitting the VCF file into replicates 1 and 2."
rep1_samples=$(bcftools query -l $purged_VCF | grep "rep1" | tr '\n' ',' | sed -e 's/,$//')
rep2_samples=$(bcftools query -l $purged_VCF | grep "rep2" | tr '\n' ',' | sed -e 's/,$//')
# Create separate VCF files for replicates 1 and 2
bcftools view -s $rep1_samples $purged_VCF > 03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.vcf.gz
bcftools view -s $rep2_samples $purged_VCF > 03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.vcf.gz

# # Tidy  up
# rm $old_header $new_header $samples_to_keep 
rm $vcf_with_all_samples $purged_VCF
