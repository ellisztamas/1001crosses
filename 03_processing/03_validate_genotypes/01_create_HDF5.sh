#!/usr/bin/env bash

# Prepare HDF5 files for QC.
#
# Filter SNPs in VCF files of the parents and F8s for those found in genes and at
# fairly high frequency. Convert to HDF5 using the Python package scikit-allel.
# 
# Input:
#    Zipped VCF files for the parents and F8s
#    GFF giving coordinates of annotated genes in TAIR10
# Output:
#    A series of intermediate text files with SNP positions
#    Zipped VCF files for the parents and F8s with fewer SNPs
#    HDF5 versions of the output VCF files for the parents and offspring
#
# Tom Ellis, 2nd August 2024

# SLURM
#SBATCH --job-name=prepare_validation_files
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=40GB
#SBATCH --qos=short
#SBATCH --time=2:00:00

# Set working directory and load the conda environment
source setup.sh

# === Input === #

# VCF files for the parents and F8s
parental_vcf=03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix_mac20.vcf.gz
progeny_vcf=03_processing/02_original_sample_sheet/output/F8_snp_matrix.vcf.gz
# VCF file for the extended RegMap set
regmap_vcf=01_data/03_parental_genotypes/Arabidopsis_2029_Maf001_Filter80.vcf.gz

# Annotated positions of genes.
genes_gff=01_data/01_reference_genome/TAIR10_GFF3_genes.gff

# === Output === #

outdir=$scratchdir/03_validate_genotypes/01_create_HDF5
mkdir -p $outdir

# Text files of SNPs in the F8s
all_snps_txt=$outdir/all_snps.txt
all_snps_bed=$outdir/all_snps.bed
intersect_snps=$outdir/snps_in_genes.tsv.gz # only those SNPs found in genes

# BED file version of the GFF file
genes_bed=$outdir/TAIR10_GFF3_genes.bed

# VCF files containing only the SNPs for QC
subset_parents=$outdir/parents_only_genic_SNPs_mac160.vcf.gz
subset_progeny=$outdir/progeny_only_genic_SNPs_mac160.vcf.gz

# HDF5 files for the parents and progeny
parents_hdf5=$outdir/parents_only_genic_SNPs_mac160.hdf5
progeny_hdf5=$outdir/progeny_only_genic_SNPs_mac160.hdf5


# === Script === #

# Prepare text files needed to subset SNPs
# List of SNPs present in the F8s, filtered for minor-allele count of least 80
bcftools view --min-ac 160 $progeny_vcf | 
    bcftools query -f"%CHROM\t%POS\n"  > $all_snps_txt
# Convert to BED. Positions are given as Chr, Position, Position +1
awk '{print $1"\t"$2"\t"$2+1}' $all_snps_txt > $all_snps_bed
# Select only features labelled as coding genes in the GFF, and convert to BED
grep "Note=protein_coding_gene" $genes_gff | awk -F'\t' '{print $1"\t"$4"\t"$5}' > $genes_bed

# Get only those SNPs that are inside genes.
# Bedtools intersect finds the SNPs inside genes, and awk reformats for bcftools
# Example output: "Chr1 6063"
bedtools intersect -a $all_snps_bed -b $genes_bed -wb | awk -F'\t' '{print $1"\t"$2}' | bgzip > $intersect_snps

# Subset the VCF files
echo "Subsetting parents"
bcftools view \
    -R $intersect_snps \
    -oZ -o $subset_parents \
    $parental_vcf
tabix -f $subset_parents

echo "Subsetting F8s"
bcftools view \
    -R $intersect_snps \
    -oZ -o $subset_progeny \
    $progeny_vcf
tabix -f $subset_progeny

# Convert the VCF files to HDF5.
echo "Converting parental VCF to HDF5."
02_library/vcf_to_HDF5.py --input $subset_parents --output $parents_hdf5

echo "Converting progeny VCF to HDF5."
02_library/vcf_to_HDF5.py --input $subset_progeny --output $progeny_hdf5


# Small VCF files for testing

# Select some interesting progeny
# 1_A3 == 9413x8423 F8 rep1 # both correct
# 1_B1 == 9413x8423 F8 rep2 # 9413 appears as 8427
# 1_E6 == 6193x6133 F8 rep1 # Pieters thinks 6133 was swapped with 6113
# 4_F11 == 6131x6113 F8 rep2 # both parents correct
# 5_A5 == 6131x6113 F8 rep1 # neither parent matches
# bcftools view \
#     --regions Chr1 \
#     --samples 1_A3,1_B1,1_E6,4_F11,5_A5 \
#     --output $outdir/test_progeny.vcf.gz \
#     $subset_progeny

# # Get the parents as well
# bcftools view \
#     --regions Chr1 \
#     --samples 6113,6131,6133,6193,8423,8427,9413 \
#     --output $outdir/test_parents.vcf.gz \
#     $subset_parents

# 02_library/vcf_to_HDF5.py --input $outdir/test_parents.vcf.gz --output $outdir/test_parents.hdf5
