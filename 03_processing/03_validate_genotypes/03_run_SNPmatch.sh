#!/usr/bin/env bash

# Compare expected genotypes of the F8s to genotypes of the parents.
#
# This takes the aggregate VCF of all F8 samples, extracts data on a single sample,
# and runs that sample through SNPmatch.
# More on SNPmatch: https://github.com/Gregor-Mendel-Institute/SNPmatch
#
# Input:
#    VCF for all F8s
#    SNPmatch database files for 1163 accessions, including all but one parent.
#    Sample sheet giving the position of each sample in the plates, expected
#        genotype, and filename.
# Output:
#    Standard SNPmatch output files.
#
# Tom Ellis, adapting code by Pieter Clauw, 1st August 2024.

# SLURM
#SBATCH --job-name=run_SNPmatch
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=3GB
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --array=291,376,342,357,299,85,86,89,81,56,28,92,20,246,211,252,267,412#0-429



# Set working directory and load SNPmatch
scratchdir=/scratch-cbe/users/$(whoami)/crosses
ml snpmatch/3.0.1-foss-2018b-python-2.7.15
ml bcftools/1.9-foss-2018b

# == Input ===
i=$SLURM_ARRAY_TASK_ID

# Directory with the SNPmatch database files
snpmatch_db=$scratchdir/03_validate_genotypes/02_SNPmatch_database
# VCF files with SNP calls for each sample
progeny_vcf=$scratchdir/03_validate_genotypes/01_create_HDF5/progeny_only_genic_SNPs_mac160.vcf.gz

# Sample sheet linking sequencing ID, and sample ID
sample_sheet=01_data/02_F8_unaligned_bams/sequencing_plates_original.csv 

# === Output === 

# Directory for the output
scratchdir=$scratchdir/03_validate_genotypes/03_run_SNPmatch/$sample_name
mkdir -p $scratchdir

# === Script === #


# Sample name for this job
sample_name=$(bcftools query -l $progeny_vcf | sed -n "${i}p")

# VCF file for a single sample
bcftools view \
    -s $sample_name \
    -o $scratchdir/$sample_name.vcf \
    $progeny_vcf
if [ $? -eq 0 ] ; then echo "Created VCF file for individual sample."; fi

# Run SNPmatch on the VCF file
snpmatch cross \
    -i $scratchdir/$sample_name.vcf \
    -d $snpmatch_db/parental_snp_matrix.hdf5 \
    -e $snpmatch_db/parental_snp_matrix.acc.hdf5 \
    -b 200000 \
    -v \
    -o $scratchdir/$sample_name
if [ $? -eq 0 ] ; then echo "SNPmatch completed successfully"; fi

# Clean up
rm $scratchdir/$sample_name.vcf

