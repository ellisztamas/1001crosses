#!/usr/bin/env bash

# Run Kruskall-Wallis GWAS seed size for North-South intercross lines.
# Plot a Manhattan plot and QQplot from the results.
#
# Tom Ellis, 28th October 2025

# SLURM
#SBATCH --job-name=08_kruskall_wallis_gwas
#SBATCH --output=05_results/23_admixed_only_GWAS/slurm/%x-%a.out
#SBATCH --error=05_results/23_admixed_only_GWAS/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=medium
#SBATCH --time=4:00:00
#SBATCH --array=4-5

# Load conda environment
source setup.sh 
set -e

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory with flowering-time BLUPs
pheno_file_array=(05_results/23_admixed_only_GWAS/output/01_subset_lines/*tsv)
phenotype_file=${pheno_file_array[$i]}
phenotype_name=${phenotype_name}

# Genotype data
genotype_hdf5=03_processing/09_impute_haplotypes/output/F8_imputed.hdf5



# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=05_results/23_admixed_only_GWAS/output/08_kruskall_wallis_gwas
mkdir -p $outdir

# Output file from GWAS
gwas_results=${outdir}/$(basename -s .tsv ${phenotype_file})_results.tsv


# === Main === #

date
echo "Running Kruskall-Wallis GWAS for phenotype file ${phenotype_file}."
echo "Using genotype data from ${genotype_hdf5}."

# Run Kruskall-Wallis GWAS
python 02_library/kruskall_wallis_GWAS.py \
    --hdf5_file $genotype_hdf5 \
    --pheno $phenotype_file \
    --output $gwas_results

# Create some plots
# Manhattan and QQ plots as .png files, which are quick
02_library/plot_gwas.py \
  --input $gwas_results \
  --outDir $outdir