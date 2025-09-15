#!/usr/bin/env bash
# 
# Identify a set of reliable SNPs in the parents.
#
# Generate a list of SNP positions for imputation:
#    - Check whether there are parents with long stretches of heterozygosity
#    - Ignoring those stretches, remove SNPs showing any heterozygosity
#    - Take only SNPs in genes, as these are largely syntenic
#
# Tom Ellis, 1st August 2025

# SLURM
#SBATCH --job-name=02_identify_reliable_SNPs
#SBATCH --output=03_processing/09_impute_haplotypes/slurm/%x.out
#SBATCH --error=03_processing/09_impute_haplotypes/slurm/%x.err
#SBATCH --mem=1GB
#SBATCH --qos=short
#SBATCH --time=8:00:00

# Set working directory and load the conda environment
source setup.sh

# === Input ===

# .ped file giving lines to be resequenced
imputation_sample_sheet=03_processing/09_impute_haplotypes/output/imputation_sample_sheet.ped
# Text file giving names of parental lines I think have large segregating chunks
segregating_parents=01_data/03_parental_genotypes/segregating_parents.txt

# Genotype data for 1163 accessions
parents_vcf=03_processing/01_parental_SNP_matrix/output/1163g.filtered_by_site.vcf.gz

# Annotated positions of genes.
genes_gff=01_data/01_reference_genome/TAIR10_GFF3_genes.gff

# === Output ===

outdir=$scratchdir/09_impute_haplotypes/02_identify_reliable_SNPs
mkdir -p $outdir

# Text file listing lines that probably are not segregating
homozygous_lines=$outdir/homozygous_lines.txt

# Text files of SNPs, removing likely pseudoheterozygous positions
homozygous_snps_txt=$outdir/homozygous_snps.txt
homozygous_snps_bed=$outdir/homozygous_snps.bed

# BED file version of the GFF file containing gene positions
genes_bed=$outdir/TAIR10_GFF3_genes.bed

# only those SNPs found in genes
intersect_snps=$outdir/homozygous_snps_in_genes.tsv.gz 


# === Main ===

# Create a file of parental lines to keep, with one name per row
# This combines the first two columns of the sample sheet, then finds and sorts
# unique elements.
# The grep command then excludes samples likely to be segregating
cut -f3,4 $imputation_sample_sheet | tr '\t' '\n' | sort -u | \
    grep -v -F -f $segregating_parents > $homozygous_lines

# Remove lines with likely-real heterozygosity, then filter for any heterozygous
# sites and minor-allele frequency
echo "Creating a list of SNPs that are probably pseudo-heterozygous."
bcftools view \
    --samples-file $homozygous_lines \
    -O u \
    $parents_vcf | \
bcftools view \
    --min-af 0.05 \
    --genotype ^het | 
bcftools query \
    -f'%CHROM\t%POS\n' \
    > $homozygous_snps_txt

echo "Subsetting SNPs that are also inside genes."
# Convert to BED. Positions are given as Chr, Position, Position +1
awk '{print $1"\t"$2"\t"$2+1}' $homozygous_snps_txt > $homozygous_snps_bed
# Select only features labelled as coding genes in the GFF, and convert to BED
grep "Note=protein_coding_gene" $genes_gff | awk -F'\t' '{print $1"\t"$4"\t"$5}' > $genes_bed
# Get only those SNPs that are inside genes.
# Bedtools intersect finds the SNPs inside genes, and awk reformats for bcftools
# Example output: "Chr1 6063" (tab separated)
bedtools intersect -a $homozygous_snps_bed -b $genes_bed -wb | \
    awk -F'\t' '{print $1"\t"$2}' | \
    bgzip \
    > $intersect_snps

