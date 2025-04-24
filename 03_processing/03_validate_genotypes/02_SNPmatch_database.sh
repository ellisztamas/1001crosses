#!/usr/bin/env bash

# Create database files needed to run SNPmatch.
# https://github.com/Gregor-Mendel-Institute/snpmatch

# Input:
#    zipped SNP matrix for the parents
# Output:
#    database files ending .csv, .hdf5, .acc.hdf5 and .csv.json

# SLURM
#SBATCH --job-name=snpmatch_db
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=20GB
#SBATCH --qos=short
#SBATCH --time=08:00:00

# Set working directory and load SNPmatch
scratchdir=/scratch-cbe/users/$(whoami)/crosses
# ml build-env/2020
ml snpmatch/3.0.1-foss-2018b-python-2.7.15
ml bcftools/1.9-foss-2018b
# === Input === #

# SNP matrix file
parental_vcf=$scratchdir/03_validate_genotypes/01_create_HDF5/parents_only_genic_SNPs_mac160.vcf.gz

# === Output === #

# Output directory for the database files
outdir=$scratchdir/03_validate_genotypes/02_SNPmatch_database
mkdir -p $outdir

output_prefix=$outdir/parental_snp_matrix

# === Script === #

# Unzip the SNP matrix
gunzip -c $parental_vcf > $outdir/unzipped_snp_matrix.vcf

# CREATE DATABASE #
snpmatch makedb \
    -i $outdir/unzipped_snp_matrix.vcf \
    -o $output_prefix