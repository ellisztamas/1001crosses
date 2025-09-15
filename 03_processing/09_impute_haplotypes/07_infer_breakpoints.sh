#!/usr/bin/env bash

# Infer the ancestry of markers in the F8s, and infer recombination breakpoints
# using an HMM.
# 
# For a single F8, this pulls out the genotypes of the parents and F8s from 
# separate VCF files, filters them to retain reliable SNPs, and passes them 
# through several R functions to infer ancestry.
# 
#
# Tom Ellis, 7th August 2025

# SLURM
#SBATCH --job-name=07_infer_breakpoints
#SBATCH --output=03_processing/09_impute_haplotypes/slurm/%x-%a.out
#SBATCH --error=03_processing/09_impute_haplotypes/slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=7,34

# Set working directory
scratchdir=/scratch-cbe/users/$(whoami)/crosses

# Load R and VCFtools.
# I could not use the conda environment here, because vcfR conflicts with things.
module load build-env/f2022
module load bcftools/1.17-gcc-12.2.0
module load r/4.2.0-foss-2021b


set -e

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

# VCF file containing all the SNPs from the published SNP matrix
reference_panel=$scratchdir/09_impute_haplotypes/03_reference_panel/parental_lines.vcf.gz
# VCF files for the progeny
# There are two, and the sample sheet has a column indicating which to use
resequenced_vcf=03_processing/04_resequencing/output/resequenced.vcf.gz
low_coverage_vcf=03_processing/03_validate_genotypes/output/F8_filtered.vcf.gz

# Text file giving SNPs to keep
homozygous_snps_in_genes=$scratchdir/09_impute_haplotypes/02_identify_reliable_SNPs/homozygous_snps_in_genes.tsv.gz


# === Output ===

outdir=$scratchdir/09_impute_haplotypes/07_infer_breakpoints
mkdir -p $outdir

# Pedigree file created by 01_phase_parents.sh
ped_file=$outdir/${progeny_name}.ped

# VCF files for a pair of parents 
cross_parents=$outdir/${progeny_name}_parents_unphased.vcf.gz
# VCF files for a single progeny
single_progeny=$outdir/${progeny_name}_unphased.vcf.gz
# Merged VCF with parents and progeny
merge_vcf=$outdir/${progeny_name}_merged.vcf.gz


# === Main ===

echo "Processing progeny ${progeny_name} with parents ${parent1_name} and ${parent2_name}"

if [[ "$offspring_panel_name" == "resequenced" ]]; then 
    echo "Taking offspring data from the resequenced data set."
    offspring_panel=$resequenced_vcf
elif [[ "$offspring_panel_name" == "lowcoverage" ]]; then
    echo "Taking offspring data from the low-coverage data set."
    offspring_panel=$low_coverage_vcf
fi



# Parents
echo "Creating VCF file for the parents"
bcftools view \
        -s $parent1_name,$parent2_name \
        $reference_panel |
    bcftools view \
        --genotype ^het \
        --output-type z \
        --output-file $cross_parents
tabix $cross_parents

# Progeny
echo "Creating a VCF file for a single progeny individual."
bcftools view \
    --regions-file $homozygous_snps_in_genes \
    --output-type u \
    $offspring_panel |
bcftools view \
    --samples $progeny_name \
    --output-type z \
    --output-file $single_progeny
tabix $single_progeny

# Merge
echo "Merging parental and offspring into a single VCF file."
bcftools merge \
    --output-type z \
    --output $merge_vcf \
    $cross_parents $single_progeny

# Run
echo "Starting the HMM"
03_processing/09_impute_haplotypes/06_run_HMM.R \
    --vcf_file $merge_vcf \
    --parent1 $parent1_name \
    --parent2 $parent2_name \
    --progeny $progeny_name \
    --mapping_error 0.01 \
    --outdir $outdir


