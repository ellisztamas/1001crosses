#!/usr/bin/env bash

# Use PLINK to summarise decay in LD for the two replicates of F8s plus the
# parents

# Input: Three VCF files
# Output: A zipped file giving LD between each SNP for each VCF file. Also three 
#     additional PLINK files.
#
# Tom Ellis, 15th December 2023

# SLURM
#SBATCH --job-name=ld_decay
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-2

date 

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

# Calculate LD for replicate 1
plink \
    --vcf ${vcf_files[$i]} \
    --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.05 --geno 0.1 --mind 0.5 \
    --thin 0.1 -r2 --ld-window 100 --ld-window-kb 100 \
    --ld-window-r2 0 \
    --out $outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )

date