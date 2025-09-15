#!/usr/bin/env bash

# Filter variants from the 1001 genomes dataset.
# 
# Filter variants based on the results of 01_data_inspection.sh and the plotting
# script.
# 
# Tom Ellis, adapting code from https://speciationgenomics.github.io/filtering_vcfs/
# 3rd January 2024

# SLURM
#SBATCH --job-name=03_filter_1001g
#SBATCH --output=03_processing/01_parental_SNP_matrix/slurm/%x.out
#SBATCH --error=03_processing/01_parental_SNP_matrix/slurm/%x.err
#SBATCH --mem=1GB
#SBATCH --qos=medium
#SBATCH --time=12:00:00

source setup.sh

# === Input ===

# Input VCF file
vcf_1163=01_data/03_parental_genotypes/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz
# List of sites to be purged
heterozygous_sites=03_processing/01_parental_SNP_matrix/summary_stats_intitial/heterozygous_sites_to_purge.tsv


# === Output ===

# Directrory to store the VCF files created
outdir=$scratchdir/01_parental_SNP_matrix/03_filter_VCF
mkdir -p $outdir

# Target for an intermediate VCF file
vcf_without_het_sites=$outdir/vcf_without_het_sites.vcf

# Target path for the output
vcf_out=03_processing/01_parental_SNP_matrix/output/1163g.filtered_by_site.vcf.gz

# Text file giving positions of filtered SNPs for later.
variable_sites=03_processing/01_parental_SNP_matrix/output/biallelic_snp_positions.tsv.gz

# === Main ===

echo "Filtering out heterozygous sites."
bcftools view -T ^$heterozygous_sites $vcf_1163 > $vcf_without_het_sites

echo "Filtering for mac, missingness, quality and depth."
missingness=0.9
min_quality=30
min_depth=10
max_depth=40


vcftools \
        --gzvcf $vcf_without_het_sites  \
        --remove-indels \
        --max-missing $missingness \
        --minQ $min_quality \
        --min-meanDP $min_depth \
        --max-meanDP $max_depth \
        --remove-indels \
        --min-alleles 2 \
        --max-alleles 2 \
        --recode \
        --recode-INFO-all \
        --stdout | \
    bgzip -c > $vcf_out
tabix $vcf_out

# create targets file
echo "Creating targets file."
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $vcf_out | bgzip -c > $variable_sites
tabix -s1 -b2 -e2 $variable_sites

cp $vcf_out 03_processing/01_parental_SNP_matrix/output/