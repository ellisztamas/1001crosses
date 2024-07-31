#!/usr/bin/env bash

# Use PLINK to create a matrix of pairwise LD for sites across the genome, 
# including between chromosomes

# Input: Three VCF files
# Output: A zipped file giving LD between each SNP for each VCF file. Also three 
#     additional PLINK files.
#
# Tom Ellis, 17th July 2024

# SLURM
#SBATCH --job-name=ld_gw
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=40GB
#SBATCH --qos=medium
#SBATCH --time=1-00:00:00
#SBATCH --array=0-2

source setup.sh

# === Input files ===

i=$SLURM_ARRAY_TASK_ID

# Path to a VCF files with dubious samples removed
rep1=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.vcf.gz
rep2=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.vcf.gz
# VCF file for the parents
parents=03_processing/04_pieters_VCF/output/parental_snp_matrix.vcf.gz
# Make a Bash array of the three VCF files
vcf_files=($rep1 $rep2 $parents)

# === Output === 

# Output directory
outdir=05_results/02_ld_decay/output
mkdir -p $outdir

output_prefix=$outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )_ldmatrix

# === Script === #

plink \
    --vcf ${vcf_files[$i]} \
    --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.05 --geno 0.1 --mind 0.5 \
    -r2 gz \
    --inter-chr \
    --ld-window-r2 0 --thin 0.1 \ 
    --out $output_prefix