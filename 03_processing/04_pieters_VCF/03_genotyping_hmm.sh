#!/usr/bin/env bash

# Use SNPmatch to call parental genotypes across an F2 genome.

# Inputs:
#     Zipped VCF files for the two replicates, excluding samples that could not
#         be validated, and with sample names as headers.
# Outputs:
#     

# Tom Ellis, 21st June 2024

# SLURM
#SBATCH --job-name=hmm_genotyping
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=8:00:00

date 

source setup.sh

# == Inputs ==

# VCF files for the two sets of F8s
rep1_vcf=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.vcf.gz
rep2_vcf=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.vcf.gz

plantID=6255x6136_rep1

snpmatch_db=/groups/nordborg/projects/epiclines/002.pedigree/03_processing/05_genotyping/output/snpmatch_database/snpmatch_db_parents_only.acc.hdf5

# === Output files ===

outdir=03_processing/04_pieters_VCF/output/hmm_genotyping
mkdir -p $outdir

# === Script ===

snpmatch genotype_cross \
    -v \
    -e $snpmatch_db \
    -p "6255x6136" \
    -i $rep1_vcf \
    -o $outdir/${plantID}.tsv \
    -b 10000