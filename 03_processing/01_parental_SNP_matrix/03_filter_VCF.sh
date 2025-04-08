#!/usr/bin/env bash

# Filter variants from the parental VCF.
# 
# Filter variants based on the results of 01_data_inspection.sh and the plotting
# script.
# 
# Input:
    # - Full VCF file from Brachi et al.
    # - Text file of heterozygous variants to purge created by 01_inspect_initial_snp_matrix.sh
# Output: A zipped VCF file with variants filtered for indels, minor allele count,
#     missing data, sequence quality and depth.
# 
# Tom Ellis, adapting code from https://speciationgenomics.github.io/filtering_vcfs/
# 3rd January 2024

# SLURM
#SBATCH --job-name=filter_parental_VCF
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=40GB
#SBATCH --qos=medium
#SBATCH --time=24:00:00

source setup.sh

# Input VCF file
vcf_in=01_data/03_parental_genotypes/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz
# List of sites to be purged
heterozygous_sites=03_processing/01_parental_SNP_matrix/summary_stats_intitial/heterozygous_sites_to_purge.tsv

# Directrory to store the VCF files created
outdir=03_processing/01_parental_SNP_matrix/output
# Target for an intermediate VCF file
vcf_intermediate=$outdir/temp_parental_SNP_matrix.vcf.gz
# Target path for the output
vcf_out=03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix_mac20.vcf.gz

# set filters
mac=20
missingness=0.9
min_quality=30
min_depth=10
max_depth=40

# Use bcftools to filter heterozygous sites.
bcftools view -T ^$heterozygous_sites $vcf_in > $vcf_intermediate
if [ $? -eq 0 ] ; then echo "Filtering heterozygous sites completed successfully"; fi

# Use vcftools to filter mac, missingness, quality an depth.
vcftools --gzvcf $vcf_intermediate  \
    --remove-indels \
    --mac $mac \
    --max-missing $missingness \
    --minQ $min_quality \
    --min-meanDP $min_depth \
    --max-meanDP $max_depth \
    --recode --stdout | \
bgzip -c > $vcf_out
if [ $? -eq 0 ] ; then echo "Filtering with vcftools completed successfully"; fi

tabix $vcf_out

rm $vcf_intermediate