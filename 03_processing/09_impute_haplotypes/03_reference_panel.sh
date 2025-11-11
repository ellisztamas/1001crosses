#!/usr/bin/env bash

# Create a VCF file for the parents only.

# Most parents are taken from the 1163 accessions in the 1001genomes dataset.
# 1435, 5835 and 6199 are taken from the RegMap panel because I think the genotypes
# in the 1001genomes dataset are incorrect.
# 1137 and 1074 are taken from the resequenced lines because they are missing 
# from either the 1001genomes or RegMap panels.

# This begins by subsetting the required accessions from the three VCF files, then
# merges these.

# Note that the RegMap panel contains phased, imputed genotypes but does not
# include genotype likelihoods (GL fields).
#
# Tom Ellis, 9th April 2025

# SLURM
#SBATCH --job-name=03_reference_panel
#SBATCH --output=03_processing/09_impute_haplotypes/slurm/%x.out
#SBATCH --error=03_processing/09_impute_haplotypes/slurm/%x.err
#SBATCH --mem=1GB
#SBATCH --qos=short
#SBATCH --time=8:00:00

# Set working directory and load the conda environment
source setup.sh


# === Input ===

# VCF files with parental data
# 1163 acessions containing most of the parents
g1163=03_processing/01_parental_SNP_matrix/output/1163g.filtered_by_site.vcf.gz
# Regmap panel, with probably-correct genotypes for 1435 and 5358
regmap=03_processing/01_parental_SNP_matrix/output/RegMap_overlapping_with_1163g.vcf.gz
# Resequenced lines, including 1137 and 1074
resequenced=03_processing/04_resequencing/output/resequenced.vcf.gz

# .ped file giving lines to be resequenced
sample_sheet=03_processing/09_impute_haplotypes/output/imputation_sample_sheet.ped

# Text file giving SNP positions to use
homozygous_snps_in_genes=$scratchdir/09_impute_haplotypes/02_identify_reliable_SNPs/homozygous_snps_in_genes.tsv.gz 

# === Output ===

outdir=$scratchdir/09_impute_haplotypes/03_reference_panel
mkdir -p $outdir

# Intermediate VCF files to be merge later
# Unique parents only, excluding 1435, 5835, 6199, 1137 and 1074
most_parents=$outdir/most_parents.vcf.gz
# Two accessions that were likely incorrect in the 1001genomes dataset
incorrect_parents=$outdir/incorrect_parents.vcf.gz
# Two accessions that were missing entirely and we sequenced de novo
missing_parents=$outdir/missing_parents.vcf.gz

# Text file giving unique parents to extract from the main VCF file.
list_of_most_parents=$outdir/list_of_most_parents.txt

# VCF file with parents merged.
merged_vcf=$outdir/parental_lines.vcf.gz
# Final VCF file with all parents
reference_panel=$outdir/parental_lines_genic_SNPs.vcf.gz
# Text file with the names of the parental lines
line_names=$outdir/parental_line_names.txt


# === Main ===

# This excludes 1435, 6199 and 5835, which we will get from the RegMap panel
echo "Extracting the unique parents from the main VCF file."
cut -f3,4 $sample_sheet | \tr '\t' '\n' | sort -u |
    grep -v "1435\|5835\|6199" \
     > $list_of_most_parents
bcftools view \
    -S $list_of_most_parents \
    -O z \
    -o $most_parents \
    $g1163
tabix $most_parents

echo "Extracting three accessions that were likely incorrect or missing in the 1001genomes dataset."
bcftools view \
    -s 1435,5835,6199 \
    -O z \
    -o $incorrect_parents \
    $regmap
tabix $incorrect_parents

echo "Extracting two accessions that were missing entirely and we sequenced de novo."
bcftools view \
    -s 1137,1074 \
    -O z \
    -o $missing_parents \
    $resequenced
tabix $missing_parents

# Merge the three VCF files into a single VCF file
# Keep only biallelic SNPs
# This will be used as a reference panel for the Beagle imputation
echo "Merging the three VCF files into a single VCF file."
bcftools merge \
    $most_parents $incorrect_parents $missing_parents |
  bcftools annotate \
    --set-id '%CHROM\_%POS' \
    -O z \
    -o $merged_vcf
tabix $merged_vcf

 03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz -Oz -o parental_lines.vcf.gz

echo "Filtering for biallelic SNPs in genes."
bcftools view \
        --regions-file $homozygous_snps_in_genes \
        --min-ac 2:minor \
        --min-alleles 2 \
        --max-alleles 2 \
        --types snps \
        --output-type z \
        --output $reference_panel \
        $merged_vcf
tabix $reference_panel

echo "Creating a text file with the names of the parental lines."
bcftools query -l $reference_panel > $line_names

echo "Copying the reference panel and line names to the output directory."
cp ${merged_vcf}* 03_processing/09_impute_haplotypes/output/
cp ${reference_panel}* 03_processing/09_impute_haplotypes/output/
cp $line_names 03_processing/09_impute_haplotypes/output/