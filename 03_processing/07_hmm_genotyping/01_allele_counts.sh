#!/usr/bin/env bash

# Prepare allele count files for RTIGER
#
# For each validated F8 genotype create a text file giving 
#
# 1 	SeqID 	String 	Chromosme ID
# 2 	Pos 	init(>=0) 	Position of the SNP marker
# 3 	RefA 	char 	Reference allele
# 4 	RefC 	int(>=0) 	Number of reads with reference alele
# 5 	AltA 	char 	Alternate allele
# 6 	AltF 	int(>=0) 	Number of reads with alternate allele
#
# Tom Ellis, 6th March 2024

# SLURM
#SBATCH --job-name=allele_counts
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=2:00:00

# Set working directory and load the conda environment
source setup.sh

# === Input ===

i=1

ibdpainting_results=03_processing/03_validate_genotypes/output/ibdpainting_results.csv

parents_vcf=03_processing/03_validate_genotypes/output/parents_only_genic_SNPs_mac160.vcf.gz
progeny_vcf=03_processing/03_validate_genotypes/output/F8_filtered.vcf.gz

# === Output ===

outdir=$workdir/07_hmm_genotyping/01_allele_counts
mkdir -p $outdir

sample_list=$outdir/sample_list.txt

rtiger_table_progeny=$outdir/rtiger_table_progeny.txt
rtiger_table_parent1=$outdir/rtiger_table_parent1.txt
rtiger_table_parent2=$outdir/rtiger_table_parent2.txt

# === Main === 

# Pull out the F8s that were validated as correct
grep -e "[0-9]x[0-9]|" -e ",correct," $ibdpainting_results | \
    awk '{gsub(/ F8 /, "_", $4); print $4}'  FS=, \
    > $sample_list

# Sample names to extract
progeny_name=sed -n "${i}p" < $sample_list 
parent1_name=$(echo $progeny_name | cut -d'x' -f1)
parent2_name=$(echo $progeny_name | cut -d'x' -f2 | cut -d'_' -f1)

# Extract a table of chr, position, ref allele, ref allele count, alt allele, alt allele count
# [\t%AD] pulls out counts to both alleles separated by a comma (e.g. 15,0)
# and appears twice.
# The awk commands subset each of these to extract only what is on the left and
# right of the comma. 
bcftools query \
    -s $progeny_name \
    -f '%CHROM\t%POS\t%REF[\t%AD]\t%ALT[\t%AD]\n' \
    $progeny_vcf |
    awk -F '\t' '{
        gsub(/,[0-9]+/, "", $4);
        gsub(/[0-9]+,/, "", $6);
        print
        }' \
    > $rtiger_table_progeny

# Get the allele count for parent 1
bcftools query \
    -s $parent1_name \
    -f '%CHROM\t%POS\t%REF[\t%AD]\t%ALT[\t%AD]\n' \
    $parents_vcf |
    awk -F '\t' '{
        gsub(/,[0-9]+/, "", $4);
        gsub(/[0-9]+,/, "", $6);
        print
        }' \
    > $rtiger_table_parent1

# Get the allele count for parent 2
bcftools query \
    -s $parent2_name \
    -f '%CHROM\t%POS\t%REF[\t%AD]\t%ALT[\t%AD]\n' \
    $parents_vcf |
    awk -F '\t' '{
        gsub(/,[0-9]+/, "", $4);
        gsub(/[0-9]+,/, "", $6);
        print
        }' \
    > $rtiger_table_parent2