#!/usr/bin/env bash

# Create HDF5 versions of the parental and F8 snp matrices.
#
# Split the VCF file with all valid F8 genotypes into cohorts 1 and 2, then 
# convert all four VCF files to HDF5.
#
# Tom Ellis, 11th April 2025

# SLURM
#SBATCH --job-name=06_split_and_hdf5
#SBATCH --output=03_processing/05_imputation/slurm/%x.out
#SBATCH --error=03_processing/05_imputation/slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=4:00:00

# Set working directory and load the conda environment
source setup.sh

# === Input ===

parents_vcf=$scratchdir/05_imputation/02_reference_panel/parental_lines.vcf.gz
progeny_vcf=$scratchdir/05_imputation/04_merge_and_hdf5/F8_phased_imputed.vcf.gz


# === Output ===

outdir=$scratchdir/05_imputation/${SLURM_JOB_NAME}
mkdir -p $outdir

# Text files with line names for the two F8 cohorts.
cohort1_names=$outdir/cohort1_names.txt
cohort2_names=$outdir/cohort2_names.txt

# VCF files splitting F8s by replicate cohort
cohort1_vcf=$outdir/F8_cohort1_phased_imputed.vcf.gz
cohort2_vcf=$outdir/F8_cohort2_phased_imputed.vcf.gz

# Hdf5 files
parents_hdf5=$outdir/$(basename -s .vcf.gz $parents_vcf).hdf5
progeny_hdf5=$outdir/$(basename -s .vcf.gz $progeny_vcf).hdf5
cohort1_hdf5=$outdir/F8_cohort1_phased_imputed.vcf.gz
cohort2_hdf5=$outdir/F8_cohort2_phased_imputed.vcf.gz


# === Main ===

# Create separate VCF files for the two progeny cohorts
# Text files giving line names.
bcftools query -l $progeny_vcf | grep "rep1" > $cohort1_names
bcftools query -l $progeny_vcf | grep "rep2" > $cohort2_names
# Split
bcftools view \
    -S $cohort1_names \
    -Oz \
    -o $cohort1_vcf \
    $progeny_vcf
tabix $cohort1_vcf

bcftools view \
    -S $cohort2_names \
    -Oz \
    -o $cohort2_vcf \
    $progeny_vcf
tabix $cohort2_vcf

Convert to Hdf5 
echo "Converting parental VCF to HDF5 format."
python 02_library/vcf_to_HDF5.py \
    --input  $parents_vcf \
    --output $parents_hdf5
echo "Converting VCF for all progeny to HDF5 format."
python 02_library/vcf_to_HDF5.py \
    --input  $progeny_vcf \
    --output $progeny_hdf5

echo "Converting VCF for F8 cohort 1 to HDF5 format."
python 02_library/vcf_to_HDF5.py \
    --input  $cohort1_vcf \
    --output $cohort2_hdf5

echo "Converting VCF for F8 cohort 2 to HDF5 format."
python 02_library/vcf_to_HDF5.py \
    --input  $cohort2_vcf \
    --output $cohort2_hdf5

cp ${outdir}/*hdf5 03_processing/05_imputation/output/
cp ${outdir}/*vcf.gz* 03_processing/05_imputation/output/