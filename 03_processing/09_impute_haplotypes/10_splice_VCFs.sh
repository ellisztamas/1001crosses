#!/usr/bin/env bash

# Splice up parental genotype calls at recombination breakpoints to create a new
# VCF file for each F8.
#
# Start with a BED file indicating start and end points of parental haplotype.
# Extract the relevant regions from VCF files of the parents, then concatenate 
# and sort by positions.
#
# For regions that are deemed to be heterozygous, or have private haplotypes
# (i.e. the seed stock for one of the parents was heterozygous, so the F8 does
# not match the parental genotype in the 1001-genomes data) this instead copies
# from the offspring VCF file.
#
# Tom Ellis, 11th August 2025

# SLURM
#SBATCH --job-name=10_splice_VCF
#SBATCH --output=03_processing/09_impute_haplotypes/slurm/%x-%a.out
#SBATCH --error=03_processing/09_impute_haplotypes/slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --array=1-39,41-60,62-411

# Set working directory
source setup.sh

# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# Table of valid F8 genotypes and their parents in .ped format
imputation_sample_sheet=03_processing/09_impute_haplotypes/output/imputation_sample_sheet.ped

# Name of the F8 genotype to process
progeny_name=$(awk -F'\t' -v row=$i 'NR==row {print $2}' $imputation_sample_sheet)
# Names of the parents
parent1_name=$(awk -F'\t' -v row=$i 'NR==row {print $3}' $imputation_sample_sheet)
parent2_name=$(awk -F'\t' -v row=$i 'NR==row {print $4}' $imputation_sample_sheet)
# Value indicating whether the genotype is from the resequenced or low-coverage datasets
offspring_panel_name=$(awk -F'\t' -v row=$i 'NR==row {print $1}' $imputation_sample_sheet )

# Bed file giving start and end positions of each haplotype.
bed_file=03_processing/09_impute_haplotypes/output/09_fill_haplotype_gaps/${progeny_name}_haplotypes.bed

# VCF files for the panels of parents and progeny
parents_vcf=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
# VCF files for the progeny
# There are two, and the sample sheet has a column indicating which to use
resequenced_vcf=03_processing/04_resequencing/output/resequenced.vcf.gz
low_coverage_vcf=03_processing/03_validate_genotypes/output/F8_filtered.vcf.gz


# === Output ===

outdir=$scratchdir/09_impute_haplotypes/10_splice_VCF
mkdir -p $outdir

# Separate .bed files for parent1, parent2 and private haplotypes.
parent1_bed=$outdir/${progeny_name}_parent1.bed
parent2_bed=$outdir/${progeny_name}_parent2.bed
private_bed=$outdir/${progeny_name}_private.bed

# A text file containing only $progeny_name used to rename the sample in VCF files.
new_header=$outdir/${progeny_name}_new_header.txt

# VCF files containing haplotypes extracted from the each donor
parent1_haplotypes=$outdir/${progeny_name}_parent1_haplotypes.vcf.gz
parent2_haplotypes=$outdir/${progeny_name}_parent2_haplotypes.vcf.gz
private_haplotypes=$outdir/${progeny_name}_private_haplotypes.vcf.gz
# A VCF file containing all haplotypes
concatenated_vcf=$outdir/${progeny_name}_concatenated_haplotypes.vcf.gz
# Sorted version of the concatenated VCF
sorted_vcf=$outdir/${progeny_name}_sorted_haplotypes.vcf.gz



# === Main ===

# If the input bash file does not exist, exit immediately.
if [[ ! -f $bed_file ]]; then
    printf ".bed file not found: ${bed_file}"
    exit
fi

printf "Processing progeny ${progeny_name} with parents ${parent1_name} and ${parent2_name}.\n"
if [[ "$offspring_panel_name" == "resequenced" ]]; then 
    printf "Taking offspring data from the resequenced data set.\n"
    private_vcf=$resequenced_vcf
elif [[ "$offspring_panel_name" == "lowcoverage" ]]; then
    printf "Taking offspring data from the low-coverage data set.\n"
    private_vcf=$low_coverage_vcf
fi

# If there are only private haplotypes, then just copy from the original VCF
# file and terminate the program.
# Otherwise, subset the parental and offspring VCF files, and merge later
if ! grep -q -E "parent1|parent2" $bed_file; then
    printf "Bed file contains only private haplotypes.\n"
    printf "The original input VCF for the progeny will be used as the output VCF file.\n"
    bcftools view \
        --samples $progeny_name \
        --output-type z \
        --output $sorted_vcf \
        $private_vcf
    tabix $sorted_vcf
    exit
else
    printf "\nThe .bed file contains a mix of haplotypes.\n"
    printf "The output VCF will be made up of the consituent parental and private haplotypes.\n"

    # Text file containing only the sample name, used for reheadering.
    printf $progeny_name > $new_header

    # String containing paths to be merged.
    # If there are maternal or paternal haplotypes to merge, these will be
    # concatenated before being passed to bcftools merge.
    files_to_concatenate=$private_haplotypes

    # Extract regions of interest from each donor.
    # Parent 1
    if grep -q -E "parent1" $bed_file; then
        printf "Extracting and sorting ${parent1_name} haplotypes.\n"
        grep "parent1" $bed_file > $parent1_bed
        bcftools view \
                --samples $parent1_name \
                --regions-file $parent1_bed \
                --output-type z \
                $parents_vcf | \
            bcftools reheader \
                --samples $new_header \
                --output $parent1_haplotypes
        tabix $parent1_haplotypes
        files_to_concatenate="$files_to_concatenate $parent1_haplotypes"
    fi
    # Parent 2
    if grep -q -E "parent2" $bed_file; then
        printf "Extracting and sorting ${parent2_name} haplotypes.\n"
        grep "parent2" $bed_file > $parent2_bed
        bcftools view \
                --samples $parent2_name \
                --regions-file $parent2_bed \
                --output-type z \
                $parents_vcf | \
            bcftools reheader \
                --samples $new_header \
                --output $parent2_haplotypes
        tabix $parent2_haplotypes
        files_to_concatenate="$files_to_concatenate $parent2_haplotypes"
    fi

    # Private haplotypes.
    # It is assumed that all bed files contain at least one private haplotypes
    # because there are gaps between the parental haplotypes.
    printf "Extracting and sorting private haplotypes.\n"
    # Subset bed file for private haplotypes
    grep -e "het" -e "private" $bed_file > $private_bed
    bcftools view \
        --samples $progeny_name \
        --regions-file $private_bed \
        --output-type z \
        --output $private_haplotypes \
        $private_vcf
    tabix $private_haplotypes

    # Concatenate VCF files
    printf "Concatenating haplotype files.\n"
    bcftools concat \
        --allow-overlaps \
        --output-type z \
        --output $concatenated_vcf \
        $files_to_concatenate

    printf "Sorting the concatenating VCF file.\n"
    bcftools sort \
        --output-type z \
        --output $sorted_vcf \
        $concatenated_vcf
    tabix $sorted_vcf
fi
