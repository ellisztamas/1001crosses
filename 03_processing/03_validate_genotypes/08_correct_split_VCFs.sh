#!/usr/bin/env bash
#
# Create VCF files for replicates 1 and 2 with correct sample names
#
# The VCF file created in 03_processing/02_original_sample_sheet gives sample
# names as their sequencing plate, and position on those plates (eg 2_G3), and
# contains mistakes. 
#
# This script changes the names to sample names (eg 6201x9399 F8 rep1), corrects
# mislabelled samples, and removes samples with no data, those that were
# self-pollinated, and those that turned out to be A. arenoa. It then splits the
# VCF file into replicates 1 and 2.
#
# Tom Ellis, 13th January 2025

# SLURM
#SBATCH --job-name=correct_split_vcf
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00

# Set a working directory and load the conda environment
source setup.sh
# We also need R, which isn't in the conda environment
module load build-env/f2022
module load r/4.2.0-foss-2021b

# == Input === #

# VCF file for the F8 individuals.
F8_snp_matrix=03_processing/02_original_sample_sheet/output/F8_snp_matrix.vcf.gz

# Table of results from IBDpainting.
ibdresults=03_processing/03_validate_genotypes/output/ibdpainting_results.csv



# === Output === # 

# Directory for the output
outdir=$scratchdir/03_validate_genotypes/08_correct_split_vcf
mkdir -p $outdir

# Text files for subsetting VCF files (created by the 07_new_sample_names.R)
vcf_sample_names_as_positions=$outdir/vcf_sample_names_as_positions.txt
vcf_sample_names_as_names=$outdir/vcf_sample_names_as_names.txt
vcf_sample_names_filtered=$outdir/vcf_sample_names_filtered.txt

# VCF file for all samples, with names rather than plate positions
F8_with_sample_names=${outdir}/F8_with_sample_names.vcf.gz
# VCF files for the two replicates, excluding self-fertilised lines, and those with no SNP data.
F8_filtered=${outdir}/F8_filtered.vcf.gz



# === Script === #

# List of sample names in position format.
bcftools query -l $F8_snp_matrix > $vcf_sample_names_as_positions

# Comma-separated list of lines to exclude
Rscript 03_processing/03_validate_genotypes/07_new_sample_names.R $outdir

# Swap the sample names
bcftools reheader \
    -s $vcf_sample_names_as_names \
    -o $F8_with_sample_names \
    $F8_snp_matrix

# Filter out bad samples
bcftools view -S $vcf_sample_names_filtered -o $F8_filtered $F8_with_sample_names
tabix $F8_filtered