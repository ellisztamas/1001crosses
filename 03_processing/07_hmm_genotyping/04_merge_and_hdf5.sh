#!/usr/bin/env bash

# Merge phased, imputed VCF files for individual F8 lines.
# Create an HDF5 file using only genic SNPs that can be used for validation.
#
# Tom Ellis, 10th April 2025

# SLURM
#SBATCH --job-name=03_merge_and_hdf5
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=4:00:00

# Set working directory and load the conda environment
source setup.sh

# === Input ===

# Directory with phased VCF files for each F8 separately
indir=$workdir/07_hmm_genotyping/03_beagle

# VCF file for the parents with common genic SNPs only previously used for validation
genic_SNPs=03_processing/03_validate_genotypes/output/snps_in_genes.tsv.gz
# reference_panel=03_processing/03_validate_genotypes/output/parents_only_genic_SNPs_mac160.vcf.gz

# === Output ===

outdir=$workdir/07_hmm_genotyping/${SLURM_JOB_NAME}
# outdir=$workdir/07_hmm_genotyping/06_merge_and_hdf5
mkdir -p $outdir

# Text file with the paths to the VCF files to merge
vcf_to_merge=$outdir/vcf_to_merge.txt

# Output files merging individual samples
merged_vcf=$outdir/F8_phased_imputed.vcf.gz

# Files with only genic SNPs to be used for validation.
validation_vcf=$outdir/F8_validation.vcf.gz
validation_hdf5=$outdir/F8_validation.hdf5


# === Main ===

# Text file listing VCF files to merge
ls $indir/*_rep?_phased.vcf.gz > $vcf_to_merge
echo "Found $(cat $vcf_to_merge | wc -l) individual VCF files to merge."


echo "Merging VCF files"
bcftools merge \
    -l $vcf_to_merge \
    -O z \
    -o $merged_vcf
tabix $merged_vcf
# Copy the merged VCF file to the project directory
cp $merged_vcf 03_processing/07_hmm_genotyping/output/


# Files for validation
echo "Subsetting the VCF to include only genic SNPs for validation."
bcftools view \
    -R $genic_SNPs \
    -Oz \
    -o $validation_vcf \
    $merged_vcf

echo "Converting the subsetted VCF to HDF5 format."
python 02_library/vcf_to_HDF5.py \
    --input  $validation_vcf \
    --output $validation_hdf5