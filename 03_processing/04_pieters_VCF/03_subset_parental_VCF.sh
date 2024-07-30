#!/usr/bin/env bash

# Subsample the VCF of parental genotypes to include only those samples that are
# a parent of an F8, and only those SNPs that vary among the F8s.
# 
# Inputs:
#     SNP matrix containing 1163 accessions.
#     Text file listing samples to include created in '01_reorder_header.R'
#     Zipped tab-separated text file containing a list of SNP positions in the
#         F8 SNP matrix, created by 01_reheader_VCF.sh
# Outputs:
#     SNP matrix containing 215 accessions and a subset number of SNPs.

# Tom Ellis, 28th June 2024

# SLURM
#SBATCH --job-name=subset_parents
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=5
#SBATCH --qos=short
#SBATCH --time=8:00:00

date 

source ./setup.sh

# == Inputs ==

# Path to the VCF file for the F8s created by Pieter
input_vcf=01_data/03_parental_genotypes/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz
parents_to_keep=$workdir/04_pieters_VCF/01_reheader_VCF/parents_to_keep.txt
# List of SNPs positions
variable_positions=$workdir/04_pieters_VCF/01_reheader_VCF/variable_positions.tsv.gz

# === Output files ===

outdir=$workdir/04_pieters_VCF/03_subset_parental_VCF
mkdir -p $outdir

vcf_215_accessions=$outdir/vcf_215_accessions.vcf.gz
outfile=03_processing/04_pieters_VCF/output/parental_snp_matrix.vcf.gz

# === Script ===

echo "Subsetting by individual."
bcftools view \
    -S $parents_to_keep \
    -oZ -o $vcf_215_accessions \
    $input_vcf
tabix $vcf_215_accessions

echo "Subsetting by site"
bcftools view \
    -R $variable_positions \
    -oZ -o $outfile \
    $vcf_215_accessions
tabix $outfile