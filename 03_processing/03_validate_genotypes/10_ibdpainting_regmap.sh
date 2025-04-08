#!/usr/bin/env bash

# Calculate IBD between F8s and the 1163 candidate parents in windows across the
# genome
#
# This runs ibdpainting on the genotype data, which breaks the genome
# into windows and compares IBD between pairs of individuals within that window.
# A cross should be a mosaic of windows that exactly match one of two parents.
#
# Tom Ellis, 24th January 2025

# SLURM
#SBATCH --job-name=regmap_ibdpainting
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=medium
#SBATCH --time=4:00:00
#SBATCH --array=31,67#1-81

# Set a working directory and load the conda environment
source setup.sh

# == Input === #

i=$SLURM_ARRAY_TASK_ID

# HDF5 file for the 1163 candidate parents
regmap_hdf5=$workdir/03_validate_genotypes/09_regmap_data/regmap_set.hdf5
# HDF5 with genotypes of the F8s
progeny_hdf5=$workdir/03_validate_genotypes/01_create_HDF5/progeny_only_genic_SNPs_mac160.hdf5

# Genotype validation results
ibdpainting_results=03_processing/03_validate_genotypes/output/ibdpainting_results.csv



# === Output === # 

# Directory for the output
outdir=03_processing/03_validate_genotypes/output/ibdpainting_regmap
mkdir -p $outdir

# Text file of ibdresults for the samples that need to be checked.
dubious_results=$outdir/dubious_samples.csv



# === Script === #

# Subset the ibdpainting results to list only the samples that need checking
grep -E "swap|unknown parent|missing haplotype" $ibdpainting_results > $dubious_results

# Line name (eg 5835x9058)
line_name=$(awk -v line="$i" 'NR==line {print $4}' FS=, $dubious_results)
# Sample name in the HDF5 file (eg 2_G3)
name_in_hdf5=$(awk -v line="$i" 'NR==line {print $1 "_" $2 $3}' FS=, $dubious_results)

# Grep the line in the sample sheet for $comma_separated_name, then pull out the cross label, e.g. '1074x1137'
expected_parents="${line_name%% *}"
expected_parents=${expected_parents/x/ } # Change the 'x' to a space, e.g. '1074 1137'

# Run the program
ibdpainting \
    --input $progeny_hdf5 \
    --reference $regmap_hdf5 \
    --window_size 500000 \
    --sample_name $name_in_hdf5 \
    --expected_match $expected_parents \
    --interactive \
    --outdir $outdir
