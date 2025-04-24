#!/usr/bin/env bash

# Merge phased, imputed VCF files for individual F8 lines.
# Convert the merged VCF file to HDF5 format.
#
# Tom Ellis, 20th March 2024

# SLURM
#SBATCH --job-name=05_ibdpainting
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=2GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=1-411

# Set working directory and load the conda environment
source setup.sh

# === Input ===

i=${SLURM_ARRAY_TASK_ID}
# i=1

# VCF and HDF5 files for progeny previously identifid as correct, for genic SNPs only
merged_vcf=$scratchdir/05_imputation/04_merge_and_hdf5/F8_validation.vcf.gz
merged_hdf5=$scratchdir/05_imputation/04_merge_and_hdf5/F8_validation.hdf5

# Reference panel of 1163 parents
reference_hdf5=03_processing/03_validate_genotypes/output/parents_only_genic_SNPs_mac160.hdf5



# === Output ===

outdir=03_processing/05_imputation/output/ibdpainting
mkdir -p $outdir



# === Main ===

sample_name=$(bcftools query -l ${merged_vcf} | sed -n "${i}p")
echo "Sample name for this job: ${sample_name}."
echo ""


# Expected parents for this sample, separated by a space
expected_parents=$(echo $sample_name | cut -d'_' -f1)
expected_parents=${expected_parents/x/ }
echo "Expected parents: ${expected_parents}."
echo ""

echo "Running IBD painting."
ibdpainting \
    --input ${merged_hdf5} \
    --reference $reference_hdf5 \
    --window_size 500000 \
    --sample_name $sample_name \
    --expected_match $expected_parents \
    --no-interactive \
    --keep_ibd_table \
    --outdir $outdir