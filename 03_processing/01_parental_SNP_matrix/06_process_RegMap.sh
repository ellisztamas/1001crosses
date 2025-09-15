#!/usr/bin/env bash

# Prepare the VCF file for the RegMap dataset.
#
# Rename the contigs to Chr1, Chr2 etc and index
# 
# Tom Ellis, 9th April 2025

# SLURM
#SBATCH --job-name=06_process_regmap
#SBATCH --output=03_processing/01_parental_SNP_matrix/slurm/%x.out
#SBATCH --error=03_processing/01_parental_SNP_matrix/slurm/%x.err
#SBATCH --mem=4GB
#SBATCH --qos=short
#SBATCH --time=4:00:00

# Set working directory and load the conda environment
source setup.sh


# === Input === #

# VCF file for the extended RegMap set
regmap_vcf=01_data/03_parental_genotypes/RegMap_panel/Arabidopsis_2029_Maf001_Filter80.vcf.gz

# VCF file for the 1163 accessions in Fernando's VCF file.
# This is used to generate SNP positions to subset the RegMap panel.
g1163=03_processing/01_parental_SNP_matrix/output/1163g.filtered_by_site.vcf.gz

# Text file giving positions of filtered SNPs for later.
targets_file=03_processing/01_parental_SNP_matrix/output/biallelic_snp_positions.tsv.gz


# === Output === #

outdir=$scratchdir/01_parental_SNP_matrix/06_process_regmap
mkdir -p $outdir

# The RegMap VCF file, but in bgzip format.
regmap_bgzipped=$outdir/regmap_bgzipped.vcf.gz

# Text file stating how contigs names should be corrected
chr_names=$outdir/chr_names.txt
# RegMap VCF with corrected filenames
vcf_with_Chr=$outdir/RegMap_with_ChrNames.vcf.gz

# A table of SNP positions in the 1163 accessions.
targets_file=$outdir/targets_files.tsv.gz
# VCF file with those SNPs.
subset_vcf=$outdir/RegMap_overlapping_with_1163g.vcf.gz


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
bcftools annotate --rename-chrs $chr_names $regmap_bgzipped | bgzip -c > $vcf_with_Chr
tabix $vcf_with_Chr

echo "Subsetting the SNPs to use only those that match the 1001 genome dataset."
bcftools view \
    --targets-file $targets_file \
    --output-type z \
    --output $subset_vcf \
    $vcf_with_Chr
tabix $subset_vcf

cp ${subset_vcf}* 03_processing/01_parental_SNP_matrix/output/