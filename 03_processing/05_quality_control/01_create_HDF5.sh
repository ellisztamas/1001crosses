#!/usr/bin/env bash

# Create HDF5 files for the parents and F8s using scikit-allel
#
# Input: VCF files for the parents and F8s
# Output: HDF5 files maintaining all the files in the VCFs
#
# Tom Ellis, adapting code from http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_hdf5, 
# 19th December 2023

# SLURM
#SBATCH --job-name=create_HDF5
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=40GB
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --array=0-1
 
date

source setup.sh

# Paths to the input VCF files to be converted
parents="03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix_mac20.vcf.gz"
F8s="03_processing/03_pieters_sample_sheet/output/F8_snp_matrix_mac20.vcf.gz"
vcf_files=($parents $F8s)
# File to convert in this job
infile=${vcf_files[$SLURM_ARRAY_TASK_ID]}

# Destination for the output
outdir=03_processing/05_quality_control/output/hdf5_files
mkdir -p $outdir
outfile=$(basename ${infile/.vcf.gz/.h5} )

# Convert the VCF file.
02_library/vcf_to_HDF5.py --input $infile --output $outdir/$outfile

date