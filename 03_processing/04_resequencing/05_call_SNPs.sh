#!/usr/bin/env bash

# Call genotypes for individual samples at known variable sites.

# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023

# SLURM
#SBATCH --job-name=05_call_SNPs
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=3GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-95

# Set working directory and load conda environment
source setup.sh

# === Input files ===

i=${SLURM_ARRAY_TASK_ID}
# i=0

# Directory containing aligned BAM files
indir=$scratchdir/04_resequencing/04_align_reads
# Select one deduplicated bam file
infile_array=($(find $indir -type f -name '*dedup.bam'))
infile=${infile_array[$i]}

# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas

# Text file listing SNP positions.
variable_sites=03_processing/01_parental_SNP_matrix/output/variable_sites.tsv.gz


# === Output files ===

# Output directory
outdir=$scratchdir/04_resequencing/05_snp_calls
mkdir -p $outdir

outfile=$outdir/$(basename $infile)
outfile=${outfile/_dedup.bam/.vcf.gz}

sample_name=${outfile/.vcf.gz/.txt}

# === Script ===

echo $(basename -s '.vcf.gz' $outfile) > $sample_name

# Get genotype likelihoods, and use them to call SNPs
echo "Calling the SNPs."
bcftools mpileup \
    --min-MQ 20 \
    -a FORMAT/DP,FORMAT/AD \
    --skip-indels \
    -f $genome \
    -Ou \
    $infile | \
bcftools call \
    -m \
    --constrain alleles \
    --targets-file $variable_sites \
    -Ou | \
bcftools norm \
    -m \
    -both \
    -f $genome \
    -Oz | \
bcftools reheader \
    -s $sample_name \
    --output $outfile

tabix -f $outfile