#!/usr/bin/env bash

# Create database files needed to run SNPmatch.
# https://github.com/Gregor-Mendel-Institute/snpmatch

# Input: zipped SNP matrix for the parents
# Output: database files ending .csv, .hdf5, .acc.hdf5 and .csv.json

# SLURM
#SBATCH --job-name=snpmatch_db
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=04:00:00

# Set working directory and load conda environment
source setup.sh

# Input
# SNP matrix file
snp_matrix=03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix_mac20.vcf.gz

# Output directory for the database files
outdir=$workdir/09_snpmatch_db
mkdir -p $outdir

# Unzip the SNP matrix
gunzip -c $snp_matrix > $outdir/unzipped_snp_matrix.vcf

# CREATE DATABASE #
snpmatch makedb \
    -i $outdir/unzipped_snp_matrix.vcf \
    -o $outdir/parental_snp_matrix

# tidy up
rm $outdir/unzipped_snp_matrix.vcf