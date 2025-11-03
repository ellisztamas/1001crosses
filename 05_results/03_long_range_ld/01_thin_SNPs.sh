#!/usr/bin/env bash

# Use PLINK to subsamples SNPs so that they are no closer than 10kb to one 
# another.
#
# Tom Ellis, 23rd April 2025

# SLURM
#SBATCH --job-name=01_thin_SNPs
#SBATCH --output=05_results/03_long_range_ld/slurm/%x-%a.out
#SBATCH --error=05_results/03_long_range_ld/slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-1

date 

source setup.sh

# === Input files ===

i=$SLURM_ARRAY_TASK_ID

# VCF files for parents and F8s
parents=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
progeny=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz
# Make a Bash array of the three VCF files
vcf_files=($parents $progeny)
# File to use in this job
infile=${vcf_files[$i]}


# === Output === 

# Output directory
outdir=$scratchdir/03_long_range_ld/01_thin_SNPs
mkdir -p $outdir

# Updated VCF file with standardised ID labels
updated_vcf=$outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )_updated.vcf.gz

# Output prefix for PLINK
# Returns Plink files
# .irem  .log  .nosex  .prune.in .prune.out
output_prefix=$outdir/$(basename -s .vcf.gz ${vcf_files[$i]} )


# === Main ===

echo "Standardising ID labels so PLink parses them consistently"
bcftools annotate --set-id '%CHROM:%POS' -O z -o $updated_vcf $infile

echo "Thinning input file ${vcf_files[$i]}"
plink \
    --vcf ${updated_vcf} \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.05 \
    --geno 0.5 \
    --mind 0.5 \
    --bp-space 10000 \
    --seed 82 \
    --write-snplist \
    --out $output_prefix

echo "Copy the pruning files to the output directory"
mkdir -p 05_results/03_long_range_ld/output/
cp ${output_prefix}.snplist 05_results/03_long_range_ld/output/

