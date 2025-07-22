#!/usr/bin/env bash

# Calculate LD between pairs of SNPs that have been pruned to be in local
# linkage equilibrium.
#
# 01_thin_SNPs.sh prepares a list of SNPs that are in local linkage in 
# the parents.
# 02_LD_between_windows.py calculates the LD between pairs of SNPs in the list.
# This script runs the LD calculation in parallel using SLURM.
#
# Tom Ellis, 24th April 2025

# SLURM
#SBATCH --job-name=03_find_SNPs_in_LD
#SBATCH --output=05_results/03_long_range_ld/slurm/%x-%a.out
#SBATCH --error=05_results/03_long_range_ld/slurm/%x-%a.err
#SBATCH --mem=3GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-111

source setup.sh

# === Input files ===

i=$SLURM_ARRAY_TASK_ID

# # VCF files for parents and F8s
parents_hdf5=03_processing/05_imputation/output/parental_lines.hdf5
progeny_hdf5=03_processing/05_imputation/output/F8_phased_imputed.hdf5

# List of SNPs to compare.
snp_list=05_results/03_long_range_ld/output/parental_lines.snplist



# === Output === 

# Output directory
outdir=05_results/03_long_range_ld/output
mkdir -p $outdir

# Output files
parents_out=$outdir/parents_snps_in_LD.csv
progeny_out=$outdir/progeny_snps_in_LD.csv


# === Main ===

# Fie headers giving the locus names and a column for r2
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo "This is the first task in the array. Creating output files."
    echo "i,j,d,dprime,r2" > $parents_out
    echo "i,j,d,dprime,r2" > $progeny_out
else
    echo "This is not the first task in the array. Waiting for the first task to create output files."
    sleep 5
fi


python 05_results/03_long_range_ld/02_LD_between_windows.py \
    --genotypes $parents_hdf5 \
    --targets $snp_list \
    --chunk_size 100 \
    --chunk_ix $i \
    --outfile $parents_out

python 05_results/03_long_range_ld/02_LD_between_windows.py \
    --genotypes $progeny_hdf5 \
    --targets $snp_list \
    --chunk_size 100 \
    --chunk_ix $i \
    --outfile $progeny_out
