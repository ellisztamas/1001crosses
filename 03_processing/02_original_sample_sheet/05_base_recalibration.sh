#!/usr/bin/env bash

# Apply GATK data preprocessing steps to the aligned BAM files.
# 
# This runs BaseRecalibrator on the aligned BAMS to identify pattrns of dubious
# covariation between bases, then runs ApplyBQSR to remove them.
# BaseRecalibrator is then run on the corrected bam file, and AnalyzeCovariates 
# plots the difference before and after recalibration.
# 
# See the GATK website for more details:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
# 
# Input:
#     Directory of aligned, deduplicated and indexed BAM files including read-
    # group information
# Output:
    # Table of recalibrtaion statistics before and after filtering
    # Filtered BAM file
    # PDF plotting the difference before and after (not run)
# 
# Tom Ellis, 13th February 2024

# SLURM
#SBATCH --job-name=base_recalibration
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-428

# Terminate the script automatically if any part of it fails (returns a non-zero exit status)
set -e

# Set the working directory and activate the conda environment
source setup.sh
# Load the modules for GATK manually
ml build-env/f2021
ml gatk/4.2.6.1-gcccore-10.2.0-java-13


# Input
# Reference genome
genome=01_data/01_reference_genome/TAIR10_chr_all.fas
# VCF fie giving known sites
known_sites=03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix_mac20.vcf.gz
# Array of input files
infiles=($workdir/04_aligned_bam/*.dedup.read_groups.bam)
infile=${infiles[$SLURM_ARRAY_TASK_ID]}
infile_bs=$(basename $infile) # file name to save for later

# Output
# Directory to hold results
outdir=$workdir/05_base_recalibration
mkdir -p $outdir
# Filenames for the output files
recal_table1=$outdir/${infile_bs/.bam/_1.recal_table}
recal_table2=$outdir/${infile_bs/.bam/_2.recal_table}
recal_bam=$outdir/${infile_bs/.bam/.recal.bam}
analyse_covariates=$outdir/${infile_bs/.bam/.pdf}


gatk BaseRecalibrator \
   -I $infile \
   -R $genome \
   --known-sites $known_sites \
   -O $recal_table1
if [ $? -eq 0 ] ; then echo "Initial base recalibration run successful." ; fi

gatk ApplyBQSR \
   --bqsr-recal-file ${recal_table1} \
   -I $infile \
   -R $genome \
   -O $recal_bam
if [ $? -eq 0 ] ; then echo "ApplyBQSR finished successfully." ; fi

gatk BaseRecalibrator \
   -I $recal_bam \
   -R $genome \
   --known-sites $known_sites \
   -O $recal_table2
if [ $? -eq 0 ] ; then echo "Second base recalibration successful." ; fi

# gatk AnalyzeCovariates \
#     -before $recal_table1 \
#     -after $recal_table2 \
#     -plots $analyse_covariates
# if [ $? - eq 0 ] ; then echo "AnalyzeCovariates completed successfuly." ; fi