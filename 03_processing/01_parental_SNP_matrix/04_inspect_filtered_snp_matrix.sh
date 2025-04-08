#!/usr/bin/env bash

# Calculate summary statistics on the filtered VCF file.
#
# Input: Filtered VCF file created by 03_processing/01_parental_SNP_matrix/03_filter_VCF.sh
# Output: Text files summarising:
    # number of variants
    # heterozygosity
    # missing data
    # sequence quality
    # sequence depth
#
# Tom Ellis, adapting code from https://speciationgenomics.github.io/filtering_vcfs/
# and https://github.com/danielecook/TIL/blob/ef86a62f1aaa5c12501bde38fc3a043b7b052763/VCF/calculate_het_rate.md
#
# 3rd Janurary 2024

# SLURM
#SBATCH --job-name=inspect_filtered_snp_matrix
#SBATCH --output=03_processing/01_parental_SNP_matrix/summary_stats_filtered/%x.out
#SBATCH --error=03_processing/01_parental_SNP_matrix/summary_stats_filtered/%x.err
#SBATCH --mem=40GB
#SBATCH --qos=short
#SBATCH --time=4:00:00

# Set working directory and load the environment
source setup.sh

# Input VCF file for the parents
vcf=03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix.vcf.gz

# Output directory, including filename prefix
outdir=03_processing/01_parental_SNP_matrix/summary_stats_filtered/summary_stats_filtered
mkdir -p $outdir
# Function to calculate heterozygosity at each locus by Daniel Cook
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

# Number of SNPs
bcftools view -H $vcf | wc -l
# Minor allele frequencies
vcftools --gzvcf $vcf --freq2 --out $outdir --max-alleles 2
# Coverage depth per individual and per site
vcftools --gzvcf $vcf --depth --out $outdir
vcftools --gzvcf $vcf --site-mean-depth --out $outdir
# Quality at each site
vcftools --gzvcf $vcf --site-quality --out $outdir
# Missing data per individual and per site
vcftools --gzvcf $vcf --missing-indv --out $outdir
vcftools --gzvcf $vcf --missing-site --out $outdir
# Heterozygosity per individual and per site
vcftools --gzvcf $vcf --het --out $outdir
heterozygosity_per_site $vcf > $outdir/heterozygosity_per_site.tsv