#!/usr/bin/env bash

# Call genotypes for individual samples at known variable sites.

# Tom Ellis, 21st May 2025

# SLURM
#SBATCH --job-name=03_call_SNPs
#SBATCH --output=03_processing/07_check_parents/slurm/%x-%a.out
#SBATCH --error=03_processing/07_check_parents/slurm/%x-%a.err
#SBATCH --mem=3GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-2

# Set working directory and load conda environment
source setup.sh

# === Input files ===

i=${SLURM_ARRAY_TASK_ID}
# i=0

# Directory containing aligned BAM files
indir=$scratchdir/07_check_parents/02_align_reads/
# Select one deduplicated bam file
infile_array=($(find $indir -type f -name '*dedup.bam'))
infile=${infile_array[$i]}

# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas

# Text file listing SNP positions.
variable_sites=03_processing/01_parental_SNP_matrix/output/variable_sites.tsv.gz


# === Output files ===

# Output directory
outdir=$scratchdir/07_check_parents/03_call_SNPs
mkdir -p $outdir

# Output VCF file
outfile=$outdir/$(basename $infile)
outfile=${outfile/_val_1.fq.dedup.bam/.vcf.gz}

# A text file containing the sample name
sample_name_file=${outfile/.vcf.gz/.txt}

# === Script ===

echo "Input file: $infile"  
echo "Output file: $outfile"

# Extract only the firs t four characters of the input file name and save to a text file
sample_name=$(basename -s .vcf.gz $outfile)
echo $sample_name > $sample_name_file


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
    -s $sample_name_file \
    --output $outfile

tabix -f $outfile