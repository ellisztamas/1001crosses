#!/usr/bin/env bash

# Run Kruskall-Wallis GWAS on flowering time.
# Plot a Manhattan plot and QQplot from the results.
#
# Tom Ellis, 20th June 2025

# SLURM
#SBATCH --job-name=01_KW_GWAS
#SBATCH --output=05_results/20_kruskall_wallis_gwas/slurm/%x-%a.out
#SBATCH --error=05_results/20_kruskall_wallis_gwas/slurm/%x-%a.err
#SBATCH --mem=20GB
#SBATCH --qos=medium
#SBATCH --time=16:00:00
#SBATCH --array=0-3

# Load conda environment
source setup.sh 
set -e

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory with flowering-time BLUPs
pheno_file_array=(03_processing/06_process_phenotypes/output/*tsv)
phenotype_file=${pheno_file_array[$i]}

# HDF5 files for parents and F8s
parents_hdf5=03_processing/05_imputation/output/parental_lines.hdf5
progeny_hdf5=03_processing/05_imputation/output/F8_phased_imputed.hdf5
# Array of paths to VCF files.
# It is really important that these are in the same order as the phenotype files!
hdf5_array=($progeny_hdf5 $parents_hdf5 $progeny_hdf5 $progeny_hdf5)
input_hdf5=${hdf5_array[$i]}



# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=05_results/20_kruskall_wallis_gwas/output
mkdir -p $outdir

# Output file from GWAS
gwas_results=${outdir}/$(basename -s .tsv ${phenotype_file})_results.tsv


# === Main === #

date
echo "Running Kruskall-Wallis GWAS for phenotype file ${phenotype_file}."
echo "Using genotype data from ${input_hdf5}."

# Run Kruskall-Wallis GWAS
python 02_library/kruskall_wallis_GWAS.py \
    --hdf5_file $input_hdf5 \
    --pheno $phenotype_file \
    --output $gwas_results

# Create some plots
# Manhattan and QQ plots as .png files, which are quick
02_library/plot_gwas.py \
  --input $gwas_results \
  --outDir $outdir