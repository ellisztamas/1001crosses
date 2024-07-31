#!/usr/bin/env bash

# Create small HDF5 files for testing
#
# Subsets the VCF file for the F8s to pull out four individuals that
# SNPmatch thinks are probably correct. Two are reciprocals of one another,
# one is another good match, and the fourth is a probable match but with
# suboptimal coverage. Also subset the parent VCF file for the parents of
# these crosses.
#
# Subset the VCF files and convert to HDF5.
#
# Input: VCF files for all the parents and F8s
# Output: HDF5 files containing a subset of individuals.
#
# 20th December 2023

# SLURM
#SBATCH --job-name=create_test_HDF5
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=40GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00

date

source setup.sh

# Paths to the VCF files for parents and F8s
# Paths to the input VCF files to be converted
parents_input_vcf="03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix.vcf.gz"
F8s_input_vcf="03_processing/03_pieters_sample_sheet/output/F8_snp_matrix.vcf.gz"

# Select three F8s to sample, plus their nominal parents.
# The first two are likely to be correct.
# 5856x9452 is probably correct, but SNPmatch flagged it as dubious sequence quality. 
F8s_to_pick="6243x6180_F8_rep1,6243x6180_F8_rep2,9413x8423_F8_rep2,5856x9452_F8_rep2"
parents_to_pick="6243,6180,9413,8423,5856,9452"

# Destination paths
outdir=03_processing/05_quality_control/output/hdf5_files

# Subset each VCF file
bcftools view -s $F8s_to_pick $F8s_input_vcf         | bgzip -c > $outdir/F8_test_data.vcf.gz
bcftools view -s $parents_to_pick $parents_input_vcf | bgzip -c > $outdir/parents_test_data.vcf.gz

# Convert the VCF files to HDF5.
02_library/vcf_to_HDF5.py \
  --input $outdir/F8_test_data.vcf.gz \
  --output $outdir/F8_test_data.h5

02_library/vcf_to_HDF5.py \
  --input $outdir/parents_test_data.vcf.gz \
  --output $outdir/parents_test_data.h5

# Tidy up intermediate files
rm $outdir/*.vcf.gz

date