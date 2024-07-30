#!/usr/bin/env bash

# SNP calling, based on a SNP matrix created by Fernando Rabanal.
# The SNP calling here is done in parallel for each chromosome.
# These get merged in 08.process_VCF.sh

# Input: aligned BAM files
# Output: 
    # Text file listing bam files to process
    # Five target files for each chromosome
    # Five zipped VCF files for each chromosome.



# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023

# SLURM
#SBATCH --job-name=call_SNPs
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=10
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --array=1-5

date

# Set working directory and load conda environment
source setup.sh
# Directory containing aligned BAM files
indir=$workdir/04_aligned_bam
# generate file with location of bamfiles (one per line)
bam_list=${indir}/bam_list.txt
# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas
# SNP matrix file
training_vcf=03_processing/01_parental_SNP_matrix/output/filtered_parental_SNP_matrix_mac20.vcf.gz

# Output directory
outdir=$workdir/05_snp_calls
mkdir -p $outdir
# chromosome specific files
chr=$SLURM_ARRAY_TASK_ID
outfile=$outdir/F8_snp_matrix_chr${chr}.vcf.gz
targets_file=$outdir/targets_file${chr}.tsv.gz

# create file with bamfile paths
ls -d ${indir}/*.bam > $bam_list

# create targets file
bcftools query -r Chr${chr} -f'%CHROM\t%POS\t%REF,%ALT\n' $training_vcf | \
    awk '{gsub(/Chr/,"",$1)} {print $1, $2, $3}' OFS='\t' | \
    bgzip -c > $targets_file
tabix -s1 -b2 -e2 $targets_file

# Get genotype likelihoods, and use them to call SNPs
bcftools mpileup --min-MQ 10 -a FORMAT/DP,FORMAT/AD --skip-indels -f $genome -r $chr -b $bam_list -Ou | \
bcftools call -m --constrain alleles --targets-file $targets_file --variants-only -Oz --output $outfile
tabix $outfile

if [ $? -eq 0 ] 
then
    echo "Genotype calling completed successfuly."
fi

date