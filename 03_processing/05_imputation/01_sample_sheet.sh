#!/usr/bin/env bash

# Create a sample sheet for imputation.
#
# For samples identified as probably correct, or with obvious label mixups that 
# could be resolved, create a tab-separated pedigree file with columns:
    # 1. Whether the genotype is resequenced or low-coverage, instead of the usual familyID
    # 2. Line name
    # 3. Father
    # 4. Mother
#
# This begins by identifying samples from the VCF file for resequenced lines
#     that look correct, have segregating parents or are obvious swaps, and
#     excluding the data on parental lines 1137 and 1074 that are included in
#     that VCF file.
# It then extracts usable lines from the low-coverage data, excluding those
#     that are already given in the resequenced data.
# Where lines appear to be swapped I used updated names.
# The first column of the ped file states 'resequenced' or 'lowcoverage',
#     which can be used in a subsequent script to indicate which VCF file to
#     look in.
# Tom Ellis, 9th April 2025

source setup.sh

# === Input ===

# IBDpainting results for resequenced lines
# If these are good, use these first 
resequencing_results=03_processing/04_resequencing/output/ibdpainting_results_resequenced.tsv
# IBDpainting results for low-coverage sequencing 
lowcoverage_results=03_processing/03_validate_genotypes/output/ibdpainting_results.csv


# === Output ===

outdir=$scratchdir/05_imputation/01_sample_sheet
mkdir -p $outdir

# Ped file for resequenced lines
resequenced_ped=$outdir/resequenced.ped
# Text file listing resequenced lines, used to exclude low-coverage samples
resequenced_names=$outdir/resequenced_names.txt
# Ped file for low-coverage lines
lowcoverage_ped=$outdir/lowcoverage.ped
# Concatenate the two ped files
beagle_sample_sheet=$outdir/beagle_sample_sheet.ped


# === Main ===


# .ped file for resequenced lines
# Tab-separated file with columns:
# 1. the word 'resequenced'
# 2. Individual ID
# 3. Father ID
# 4. Mother ID
grep -e "correct\|swap\|missing haplotype" $resequencing_results | \
awk -F'\t' -v fam='resequenced' '{
        split($5,a,"_");
        split(a[1],b,"x");
        printf "%s\t%s\t%s\t%s\n", fam, $5, b[1], b[2]
    }' |
grep -v "1137\|1074" > $resequenced_ped

# Extract the names of the resequenced lines to exclude them from the low-coverage .ped file
awk -F'\t' '{print $2}' $resequenced_ped > $resequenced_names

# .ped file for low coverage lines, excluding those in the previous .ped file
# Tab-separated file with columns:
# 1. the word 'lowcoverage'
# 2. Individual ID
# 3. Father ID
# 4. Mother ID
grep -e "[0-9]x[0-9]|" -e ",correct,\|,swap,\|missing haplotype" $lowcoverage_results | \
awk -F',' -v fam='lowcoverage' '{
        split($5,a,"[ ]");         # split field 4 at spaces to separate cross from F8 and rep
        split(a[1],b,"x");         # split the cross at "x"
        printf "%s\t%s_%s\t%s\t%s\n", fam, a[1], a[3], b[1], b[2]
    }' |
grep -vf  $resequenced_names > $lowcoverage_ped

# Concatenate the two .ped files
cat $resequenced_ped $lowcoverage_ped > $beagle_sample_sheet

# Check for duplicates.
# Prints to stout
awk -F'\t' '{count[$2]++} END {for (item in count) print item "\t" count[item]}' $beagle_sample_sheet

cp $beagle_sample_sheet 03_processing/05_imputation/output/
