#!/usr/bin/env bash

# Calculate IBD between F8s and the 1163 candidate parents in windows across the
# genome
#
# This runs ibdpainting on the genotype data, which breaks the genome
# into windows and compares IBD between pairs of individuals within that window.
# A cross should be a mosaic of windows that exactly match one of two parents.
#
# Input:
#    VCF for all F8s, filtered for common SNPs in genes
#    HDF5 file for the parents created with skit-allel (see 01_create_HDF5.sh)
#        with the same SNPs as the VCF file.
#    Sample sheet giving the position of each sample in the plates, expected
#        genotype, and filename.
# Output:
#    Output files of ibdpainting
#
# Tom Ellis, 6th August 2024

# SLURM
#SBATCH --job-name=ibdpainting
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=1-429

# Set a working directory and load the conda environment
source setup.sh

# == Input === #

i=$SLURM_ARRAY_TASK_ID

# HDF5 file for the 1163 candidate parents
ref_1163=$scratchdir/03_validate_genotypes/01_create_HDF5/parents_only_genic_SNPs_mac160.hdf5
# Files with SNP calls for each sample
# There is a VCF file and an Hdf5 file with this, name and we need both, so leave
# the suffix off here.
progeny_geno=$scratchdir/03_validate_genotypes/01_create_HDF5/progeny_only_genic_SNPs_mac160

# Sample sheet linking sequencing ID, and sample ID
sample_sheet=01_data/02_F8_unaligned_bams/sequencing_plates_original.csv 



# === Output === # 

# Directory for the output
outdir=03_processing/03_validate_genotypes/output/ibdpainting
mkdir -p $outdir



# === Script === #

# Sample name for this job
sample_name=$(bcftools query -l ${progeny_geno}.vcf.gz | sed -n "${i}p")

# This requires giving the expected parents (e.g. 5835x9058), which are in $sample_sheet
# First we need to change sample names like "2_G3" to "2,G,3" 
# Extract the parts before and after the underscore
plate="${sample_name%%_*}"  # Extracts "2"
well="${sample_name#*_}"   # Extracts "G3"
comma_separated_name="${plate},${well:0:1},${well:1}" # "2,G,3"
# Grep the line in the sample sheet for $comma_separated_name, then pull out the cross label, e.g. '1074x1137'
expected_parents=$(grep $comma_separated_name $sample_sheet | grep -oE "[0-9]+x[0-9]+")
expected_parents=${expected_parents/x/ } # Change the 'x' to a space, e.g. '1074 1137'

# Run the program
ibdpainting \
    --input ${progeny_geno}.hdf5 \
    --reference $ref_1163 \
    --window_size 500000 \
    --sample_name $sample_name \
    --expected_match $expected_parents \
    --interactive \
    --outdir $outdir
