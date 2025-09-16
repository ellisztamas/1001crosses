#!/usr/bin/env bash

# Use PLINK to summarise decay in LD within 100kb
#
# Tom Ellis, 15th December 2023

# SLURM
#SBATCH --job-name=ld_decay
#SBATCH --output=05_results/02_ld_decay/slurm/%x-%a.out
#SBATCH --error=05_results/02_ld_decay/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=1

date 

source setup.sh

# === Input files ===

i=$SLURM_ARRAY_TASK_ID

# VCF files for parents and F8s
parents=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
progeny=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz
unphased=03_processing/02_original_sample_sheet/output/F8_snp_matrix.vcf.gz
strict=03_processing/08_strict_imputation/output/F8_imputed_on_genes_beagle.vcf.gz
# Make a Bash array of the three VCF files
vcf_files=($parents $progeny $unphased $strict)
# File to use in this job
infile=${vcf_files[$i]}


# === Output === 

# Output directory
outdir=05_results/02_ld_decay/output
mkdir -p $outdir

# Plink text files giving r2 and D between pairs of loci
prefix=$(basename -s .vcf.gz ${vcf_files[$i]} )
r2_file=$outdir/${prefix}_r2
d_file=$outdir/${prefix}_D


# === Main === 

# Calculate r2
plink \
    --vcf ${vcf_files[$i]} \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.05 \
    --thin 0.1 \
    --geno 0.5 \
    --r2 gz \
    --ld-window 100 \
    --ld-window-kb 100 \
    --ld-window-r2 0 \
    --out $r2_file

# Calculate D
plink \
    --vcf ${vcf_files[$i]} \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.05 \
    --thin 0.1 \
    --geno 0.5 \
    --r gz d \
    --ld-window 100 \
    --ld-window-kb 100 \
    --ld-window-r2 0 \
    --out $d_file

