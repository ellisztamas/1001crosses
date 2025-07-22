#!/usr/bin/env bash

# Permute genotypes at each locus for the parents.
#
# This applies the python script 01_permute_VCF.py to the VCF file for the 
# parents to randomise genotypes at each marker.
# The resulting files are then converted to Plink format
#
# This is done 200 times as an array job.
#
# Tom Ellis, 21st July 2025

# SLURM
#SBATCH --job-name=02_permute_parents
#SBATCH --output=05_results/22_simulated_datasets/slurm/%x-%a.out
#SBATCH --error=05_results/22_simulated_datasets/slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=1-200

# Set working directory and load the conda environment
source setup.sh


# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# VCF file to be permuted
in_vcf='03_processing/05_imputation/output/parental_lines.vcf.gz'


# === Output ===

outdir=$scratchdir/22_simulated_datasets/${SLURM_JOB_NAME}
mkdir -p $outdir

# Prefix for output files without any file extension
output_prefix=$outdir/permutd_$(printf "%03d\n" $i)

# Output VCF file with genotypes at each locus in a random order.
out_vcf=${output_prefix}.vcf


# === Main ===

# Permute genotypes
echo "Permuting genotypes for replicate ${i}."
python 05_results/22_simulated_datasets/01_permute_VCF.py \
    --input $in_vcf \
    --output $out_vcf \
    --seed $i

# Create Plink files
echo "Converting to Plink format."
plink2 \
  --vcf $out_vcf \
  --make-bed \
  --set-all-var-ids @_# \
  --out $output_prefix