#!/usr/bin/env bash

# Prepare the VCF file for the RegMap dataset.
#
# Rename the contigs to Chr1, Chr2 etc and index
# 
# Tom Ellis, 9th April 2025

# SLURM
#SBATCH --job-name=06_process_regmap
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


# === Output === #

outdir=$scratchdir/01_parental_SNP_matrix
mkdir -p $outdir

# Text file stating how contigs names should be corrected
chr_names=$outdir/chr_names.txt

regmap_bgzipped=$outdir/regmap_bgzipped.vcf.gz
# VCF with corrected filenames
vcf_with_Chr=$outdir/RegMap_with_ChrNames.vcf.gz


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

cp ${vcf_with_Chr}* 03_processing/01_parental_SNP_matrix/output/