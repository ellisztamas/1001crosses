#!/usr/bin/env bash

# SNP calling, based on a SNP matrix created by Fernando Rabanal.
# The SNP calling here is done in parallel for each chromosome.
# These get merged in 08.process_VCF.sh

# Input:
#    Directory of aligned BAM files, with corresponding .bai files
#    Reference genome
#    Matrix of SNPs in the parents
# Output:
#    'Targets file' giving SNP positions in the parents
#    five zipped VCF files.

# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023

# SLURM
#SBATCH --job-name=call_SNPs
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=medium
#SBATCH --time=1-00:00:00
#SBATCH --array=1-5

# Set working directory and load conda environment
source setup.sh

# === Input files ===

# Directory containing aligned BAM files
indir=$scratchdir/04_aligned_bam
# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas
# SNP matrix file
parental_SNP_matrix=03_processing/01_parental_SNP_matrix/output/1163g.filtered_by_site.vcf.gz
# Which chromosome to work with
chr=Chr${SLURM_ARRAY_TASK_ID}


# === Output files ===

# Output directory
outdir=$scratchdir/05_snp_calls
mkdir -p $outdir
# file with location of bamfiles (one per line)
bam_list=${outdir}/bam_list.txt
# chromosome specific files
outfile=$outdir/F8_snp_matrix_${chr}.vcf.gz
targets_file=$outdir/targets_file_${chr}.tsv.gz


# === Script ===

echo "Creating file with paths to BAM files."
find $indir/*read_groups.bam > $bam_list

# create targets file
echo "Creating targets file."
bcftools query -r ${chr} -f'%CHROM\t%POS\t%REF,%ALT\n' $parental_SNP_matrix | bgzip -c > $targets_file
tabix -s1 -b2 -e2 $targets_file

# Get genotype likelihoods, and use them to call SNPs
echo "Calling the SNPs."
bcftools mpileup --min-MQ 20 -a FORMAT/DP,FORMAT/AD --skip-indels -f $genome -r $chr -b $bam_list -Ou | \
    bcftools call -m --constrain alleles --targets-file $targets_file --variants-only -Ou | \
    bcftools norm -m -both -f $genome -Oz --output $outfile
tabix $outfile