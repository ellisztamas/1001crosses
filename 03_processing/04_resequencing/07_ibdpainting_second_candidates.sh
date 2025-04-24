#!/usr/bin/env bash

# For candidates with plausible other candidate parents, repeat ibdpainting
# plotting with these as the expected parents.
#
# Tom Ellis, 4th April 2025

# SLURM
#SBATCH --job-name=07_ibdpainting_second_candidates
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=30:00
#SBATCH --array=15

# Set working directory and load the conda environment
source setup.sh

# === Input === #

i=${SLURM_ARRAY_TASK_ID}

# Sample sheet giving samples to analyse and the candidate parents to highlight.
sample_sheet=03_processing/04_resequencing/second_candidates.tsv

# Directory containing VCF files for each sample with genic SNPs only.
# This is created by the previous script.
indir=$scratchdir/04_resequencing/06_ibdpainting/

# VCF files for the parents and F8s
ref_panel=03_processing/03_validate_genotypes/output/regmap_set.hdf5


# === Output ===

projdir=03_processing/04_resequencing/output/ibdpainting_second_candidates
mkdir -p $projdir


# === Main ===

sample_name=$(awk -F'\t' -v row=$i 'NR==row {print $4}' $sample_sheet)
expected_parents=$(awk -F'\t' -v row=$i 'NR==row {print $5 " " $6}' $sample_sheet)

echo "Validating the VCF file."
# Run the program
ibdpainting \
    --input $indir/${sample_name}_subset.hdf5 \
    --reference $ref_panel \
    --window_size 500000 \
    --sample_name $sample_name \
    --expected_match $expected_parents \
    --interactive \
    --outdir $projdir
