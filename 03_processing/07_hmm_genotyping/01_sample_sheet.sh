#!/usr/bin/env bash

# Create a sample sheet for imputation.
#
# For samples identified as probably correct, or with obvious label mixups that 
# could be resolved, create a tab-separated pedigree file with columns:
    # 1. Family ID
    # 2. Individual ID
    # 3. Father ID
    # 4. Mother

# === Input ===

# Table giving genotype validation results of IBD painting
ibdpainting_results=03_processing/03_validate_genotypes/output/ibdpainting_results.csv


# === Output ===

# Table of valid F8 genotypes and their parents in .ped format
sample_sheet=03_processing/07_hmm_genotyping/output/sample_sheet.ped



# === Main ===

# Pull out the F8s that were validated as correct
grep -e "[0-9]x[0-9]|" -e ",correct,\|,swap," $ibdpainting_results | \
awk -F',' '{
        split($5,a,"[ ]");         # split field 4 at spaces to separate cross from F8 and rep
        split(a[1],b,"x");         # split the cross at "x"
        printf "%d\t%s_%s\t%s\t%s\n", NR, a[1], a[3], b[1], b[2]
    }' > $sample_sheet