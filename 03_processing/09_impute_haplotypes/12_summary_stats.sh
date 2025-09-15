#!/usr/bin/env bash

# Calculate summary statistics about the VCF files for the imputed F8 data, and
# the observed parental data.
#
# Tom Ellis, 9th September 2025

# SLURM
#SBATCH --job-name=12_summary_stats
#SBATCH --output=03_processing/09_impute_haplotypes/slurm/%x-%a.out
#SBATCH --error=03_processing/09_impute_haplotypes/slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=0-1

# Set working directory
source setup.sh
# set -e

# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# VCF files for the offspring and parents
parents_vcf=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
progeny_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz

# Which file to use
vcf_array=($parents_vcf $progeny_vcf)
infile=${vcf_array[$i]}


# === Output ===

outdir=03_processing/09_impute_haplotypes/output/12_summary_stats
mkdir -p $outdir

# Prefix for VCF summary statistics files is the basename of the VCF file
output_prefix=$outdir/$(basename -s .vcf.gz $infile)


# === Main ===

echo "Calculating number of SNPs"
bcftools view -H $infile | wc -l

echo "Calculating minor allele frequencies."
vcftools --gzvcf $infile --freq2 --out $output_prefix --max-alleles 2

echo "Calculating coverage depth per individual"
vcftools --gzvcf $infile --depth --out $output_prefix

echo "Calculating coverage depth per site"
vcftools --gzvcf $infile --site-mean-depth --out $output_prefix

echo "Calculating missing data per individual."
vcftools --gzvcf $infile --missing-indv --out $output_prefix

echo "Calculating missing data per locus."
vcftools --gzvcf $infile --missing-site --out $output_prefix

echo "Calculating heterozygosity per individual."
vcftools --gzvcf $infile --het --out $output_prefix


Function to calculate heterozygosity at each locus by Daniel Cook
function heterozygosity_per_site {
    vcf=${1}
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[\t%SAMPLE=%GT]\n' ${vcf} | \
    awk -v OFS='\t' 'BEGIN {print "CHROM\tPOS\tREF\tALT\tn_het\tn\thet_rate"} 
    {
    n_het=gsub(/0\|1|1\|0|0\/1|1\/0|0\|2|2\|0|0\/2|2\/0|0\|3|3\|0|0\/3|3\/0|0\|4|4\|0|0\/4|4\/0/, "" ); 
    print $1, $2, $3, $4, 
    n_het,
    (NF-5),
    n_het/(NF - 5)
    }'
}
heterozygosity_per_site $infile > ${output_prefix}.lhet
