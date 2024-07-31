#!/usr/bin/env bash

# Merge VCF files of genotypes imputed by alphaPlantImpute2 for each F8 sample
# separately into a single VCF file.
# 
# Inputs:
#     400 or so VCF files, each containing a single sample.
# Outputs:
#     A single VCF file with all 400 samples.
#
# Tom Ellis, 4th July 2024
# Adapting code from https://www.biocomputix.com/post/how-to-combine-merge-vcf-bcf-files-using-bcftools-merge

# SLURM
#SBATCH --job-name=merge_imputed_genotypes
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=8:00:00

date 
set -e

source ./setup.sh

# == Inputs ==

indir=$workdir/07_hmm_genotyping/04_alphaPlantImpute

# The conversion to VCF requires a .map file, and this is as good as any
map_file=$indir/6255x6136_rep1/offspring_genotypes.map


# === Output files ===

outdir=$workdir/07_hmm_genotyping/05_merge_imputed_genotypes
mkdir -p $outdir

merged_ped=$outdir/merged_imputed_genotypes_rep1
merged_map=$outdir/merged_imputed_genotypes_rep1.map

# # Text files to store genotypes in replicates 1 and 2.
# rep1_sample_list=$outdir/rep1_sample_list.txt
# rep2_sample_list=$outdir/rep2_sample_list.txt

# Merged VCF files
rep1_vcf=$outdir/F8_rep1_SNP_matrix
rep2_vcf=$outdir/F8_rep2_SNP_matrix

# === Script ===

cat $indir/*_rep1/imputed.ped > $merged_ped

# .map file needs to be in the same directory as the .ped file, with the same prefix
cp $map_file $merged_map
# Convert to (zipped) VCF
plink2 --pedmap $merged_ped --recode vcf bgz --out $rep1_vcf

# # List of imputed VCF files for replicates 1 and 2
# find $indir/**/*_rep1.vcf.gz > $rep1_sample_list
# find $indir/**/*_rep2.vcf.gz > $rep2_sample_list

# # Merge the files
# bcftools merge --file-list $rep1_sample_list -Oz -o $rep1_vcf
# bcftools merge --file-list $rep2_sample_list -Oz -o $rep2_vcf

# tabix $rep1_vcf
# tabix $rep2_vcf