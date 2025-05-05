#!/usr/bin/env bash

# Use PLINK to prune SNPs in local LD with one another.

# This uses a greedy algorithm to go through SNPs within a window and remove 
# identify pairs that are in LD>0.5. If they are, keep the one with the higher
# MAF. This is done in a sliding window of 1,000 SNPs, with a step size of 100.
# 
# The key is that this is done on a local scale, so I can still calculate 
# long-range LD later.
#
# Tom Ellis, 23rd April 2025

# SLURM
#SBATCH --job-name=01_sequential_pruning
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0

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
outdir=$scratchdir/03_long_range_ld/${SLURM_JOB_NAME}
mkdir -p $outdir

updated_vcf=$outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )_updated.vcf.gz

output_prefix=$outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )

targets_file=$outdir/targets_file.txt

# === Main ===

# Standardise ID labels so PLink parses them consistently
bcftools annotate --set-id '%CHROM:%POS' -O z -o $updated_vcf $infile

echo "Running sequential pruning on ${vcf_files[$i]}"
# Returns Plink files
# .irem  .log  .nosex  .prune.in .prune.out
plink \
    --vcf ${updated_vcf} \
    --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.05 --geno 0.1 --mind 0.5 \
    --indep-pairwise 1000 100 0.5 \
    --out $output_prefix

awk '{print $1,"\t",$2}' FS=':' ${output_prefix}.prune.in > $targets_file

cp ${output_prefix}.prune.in 05_results/03_long_range_ld/output