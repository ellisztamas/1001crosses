#!/usr/bin/env bash

# Prepare the VCF file for the RegMap dataset.
#
# Rename the contigs, index, and create an HDF5 file
# 
# Tom Ellis, 24th January 2025

# SLURM
#SBATCH --job-name=regmap_hdf5
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=40GB
#SBATCH --qos=short
#SBATCH --time=2:00:00

# Set working directory and load the conda environment
source setup.sh

# === Input === #

# VCF file for the extended RegMap set
regmap_vcf=01_data/03_parental_genotypes/Arabidopsis_2029_Maf001_Filter80.vcf.gz

# Table of SNP positions in genes
intersect_snps=03_processing/03_validate_genotypes/output/snps_in_genes.tsv.gz

# === Output === #

outdir=$workdir/03_validate_genotypes/09_regmap_data
mkdir -p $outdir

# Text file stating how contigs names should be corrected
chr_names=$outdir/chr_names.txt

regmap_bgzipped=$outdir/regmap_bgzipped.vcf.gz
# VCF with corrected filenames
vcf_with_Chr=$outdir/regmap_set.vcf.gz
# VCF file with genic SNPs only
subset_vcf=$outdir/regmap_set_subset.vcf.gz

# HDF5 version of the VCF
regmap_hdf5=03_processing/03_validate_genotypes/output/regmap_set.hdf5


# === Main === #

# Convert to bgzipped format and index
bcftools view $regmap_vcf | bgzip -c > $regmap_bgzipped
tabix $regmap_bgzipped

# Swap chromosome labels from 1,2,3,4,5 to Chr1, Chr2, etc
# Create a space-delimited file for renaming chromosome labels in the VCF file.
echo "Renaming chromosome annotations."
cat > $chr_names << EOL
1 Chr1
2 Chr2
3 Chr3
4 Chr4
5 Chr5
EOL

echo "Renaming chromosomes."
bcftools annotate --rename-chrs $chr_names $regmap_bgzipped | 
bgzip -c > $vcf_with_Chr
tabix $vcf_with_Chr

echo "Subsetting the VCF file to include only SNPs in genes."
bcftools view \
    -R $intersect_snps \
    -oZ \
    -o $subset_vcf \
    $vcf_with_Chr
tabix $subset_vcf

echo "Converting RegMap VCF to HDF5."
02_library/vcf_to_HDF5.py --input $subset_vcf --output $regmap_hdf5
