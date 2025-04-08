#!/usr/bin/env bash

# Validate genotypes using ibdpainting

# Subset SNPs to include only genic SNPs, convert to HDF5, and validate the 
# genotypes using ibdpainting.
#
# Tom Ellis, 1st April 2025

# SLURM
#SBATCH --job-name=06_ibdpainting
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-95

# Set working directory and load the conda environment
source setup.sh

# === Input === #

i=${SLURM_ARRAY_TASK_ID}

# Input VCF file to be converted to HDF5 and validated
indir=$workdir/08_resequencing/05_snp_calls/
infile_array=($(find $indir -type f -name '*vcf.gz'))
infile=${infile_array[$i]}

# Table of SNP positions in genes
intersect_snps=03_processing/03_validate_genotypes/output/snps_in_genes.tsv.gz

# VCF files for the parents and F8s
# ref_panel=03_processing/03_validate_genotypes/output/parents_only_genic_SNPs_mac160.hdf5
ref_panel=03_processing/03_validate_genotypes/output/regmap_set.hdf5

# === Output === #

# Working directory on scratch-cbe for intermediate files
scratchdir=$workdir/08_resequencing/06_ibdpainting
mkdir -p $scratchdir

# Project directory to store IBDpainting results
projdir=03_processing/08_resequencing/output/ibdpainting
mkdir -p $projdir

# VCF file for a single sample containing genic SNPs only
subset_vcf=$scratchdir/$(basename $infile)
subset_vcf=${subset_vcf/.vcf.gz/_subset.vcf.gz}
subset_hdf5=${subset_vcf/.vcf.gz/.hdf5}


# === Main === #

echo "VCF file to be validated: ${infile}"
echo ""

# Sample and parent names
sample_name=$(basename -s '.vcf.gz' $infile)
# Grep the line in the sample sheet for $comma_separated_name, then pull out the cross label, e.g. '1074x1137'
expected_parents="${sample_name%%_*}"
expected_parents=${expected_parents/x/ } # Change the 'x' to a space, e.g. '1074 1137'

echo "Subsetting the VCF file to include only SNPs in genes."
bcftools view \
    -R $intersect_snps \
    -oZ -o $subset_vcf \
    $infile
tabix -f $subset_vcf

echo "Converting the VCF file to HDF5."
python 02_library/vcf_to_HDF5.py \
    --input $subset_vcf \
    --output $subset_hdf5

echo "Validating the VCF file."
# Run the program
ibdpainting \
    --input $subset_hdf5 \
    --reference $ref_panel \
    --window_size 500000 \
    --sample_name $sample_name \
    --expected_match $expected_parents \
    --interactive \
    --outdir $projdir