#!/usr/bin/env bash

# Use PLINK to summarise decay in LD for the two replicates of F8s plus the
# parents

# Input:
#    Three VCF files
# Output:
#    A zipped file giving LD between each SNP for each VCF file. Also three 
#        additional PLINK files.
#
# Tom Ellis, 15th December 2023

# SLURM
#SBATCH --job-name=ld_decay
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-1

date 

source setup.sh

# === Input files ===

i=$SLURM_ARRAY_TASK_ID

# VCF files for parents and F8s
parents=03_processing/05_imputation/output/parental_lines.vcf.gz
progeny=03_processing/05_imputation/output/F8_phased_imputed.vcf.gz
# Make a Bash array of the three VCF files
vcf_files=($parents $progeny)
# File to use in this job
infile=${vcf_files[$i]}


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
    --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 100 \
    --ld-window-r2 0 \
    --out $outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )

date