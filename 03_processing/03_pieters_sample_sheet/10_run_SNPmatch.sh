#!/usr/bin/env bash

# Compare expected genotypes of the F8s to genotypes of the parents.
#
# This takes the aggregate VCF of all F8 samples, extracts data on a single sample,
# and runs that sample through SNPmatch.
# More on SNPmatch: https://github.com/Gregor-Mendel-Institute/SNPmatch
#
# Input: VCF for all F8s, SNPmatch database files.
# Output: For each F8 separately, a human-readable CSV file sumamrising the JSON
#   output of SNPmatch.
#
# Tom Ellis, adapting code by Pieter Clauw, 27th November 2023.

# SLURM
#SBATCH --job-name=run_SNPmatch
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=3GB
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --array=1-429

source setup.sh

# Directory with the SNPmatch database files
snpmatch=$workdir/06_snpmatch/db
# VCF files with SNP calls for each sample
F8_snp_calls=$workdir/05_snp_calls/F8_snp_matrix_mac20.vcf.gz
# Sample name for this job
i=$SLURM_ARRAY_TASK_ID
sample_name=$(bcftools query -l $F8_snp_calls | sed -n "${i}p")
# Directory for the output
results=$workdir/06_snpmatch/$sample_name
mkdir -p $results
# Directory to stage out the results
outdir=03_processing/03_pieters_sample_sheet/output/snpmatch/
mkdir -p $outdir

# VCF file for a single sample
bcftools view \
    -s $sample_name \
    -o $results/$sample_name.vcf \
    $F8_snp_calls
if [ $? -eq 0 ] ; then echo "Created VCF file for individual sample."; fi

# Run SNPmatch on the VCF file
snpmatch cross \
    -i $results/$sample_name.vcf \
    -d $snpmatch/parental_snp_matrix.hdf5 \
    -e $snpmatch/parental_snp_matrix.acc.hdf5 \
    -b 200000 \
    -v \
    -o $results/$sample_name
if [ $? -eq 0 ] ; then echo "SNPmatch completed successfully"; fi

# Clean up
rm $results/$sample_name.vcf

# Create a human-readable CSV file from the JSON output of SNPmatch
parents=$(echo $sample_name | cut -d '_' -f 1)
python 02_library/CSV_from_SNPmatch.py \
    -i $results \
    -d $snpmatch/parental_snp_matrix.hdf5 \
    -e $parents \
    -f $sample_name \
    -o $outdir/$sample_name
if [ $? -eq 0 ] ; then echo "CSV_from_SNPmatch completed successfully."; fi