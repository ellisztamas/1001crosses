#!/usr/bin/env bash

# Pull out the genotypes for flowering time SNPs.
# 
# This takes a file of candidate SNPs and pulls out genotypes at each.
# Note that the output has SNPs as rows and samples as columns, and you probably
# want the transpose of that. Heterozygotes are converted to NA.
#
# Tom Ellis, 1st July 2025

# SLURM
#SBATCH --job-name=06_top_SNP_genotypes
#SBATCH --output=05_results/16_gemma_ft12/slurm/%x-%a.out
#SBATCH --error=05_results/16_gemma_ft12/slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-3

# Set working directory and load the conda environment
source setup.sh


# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Directory with flowering-time BLUPs
pheno_file_array=(03_processing/06_process_phenotypes/output/flowering_time*tsv)
phenotype_file=${pheno_file_array[$i]}
pheno_basename=$(basename -s .tsv $phenotype_file)

# VCF files for parents and F8s
progeny_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz
parents_vcf=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
# Array of paths to VCF files.
# It is really important that these are in the same order as the phenotype files!
vcf_array=($progeny_vcf $parents_vcf $progeny_vcf $progeny_vcf)
input_vcf=${vcf_array[$i]}

# A hand-made list of markers that look like the tips of peaks
snp_list_csv=05_results/16_gemma_ft12/output/top_SNPs/candidate_peak_positions.csv


# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=$scratchdir/16_gemma_ft12/06_top_SNP_genotypes
mkdir -p $outdir

# The list of SNPs, but with tabs as the delimiter, for BCFtools
snp_list_tsv=$outdir/candidate_peak_positions.tsv
# Text file listing line names in the phenotype file
samples_file=$outdir/${pheno_basename}_sample_file.txt

# Output files from VCF tools.
# This has a row for each marker and a column for each individual
# Genotypes are in VCF format (0|0, 0|1, ./, etc)
genotypes_GT_format=$outdir/${pheno_basename}_GT_calls.txt
# The same data, but with genotypes converted to integers (or NA)
genotypes_integer_format=${genotypes_GT_format/_GT_calls.txt/top_SNPs.tsv}


# === Main ===

# Convert SNP list to be tab-delimited
awk -F',' '{print "Chr"$1, $2}' OFS='\t' $snp_list_csv > $snp_list_tsv

# Text file with sample names
awk -F'\t' '{print $2}' $phenotype_file > $samples_file

# Pull genotypes at focal loci in VCF format
echo "Creating a matrix if GT calls"
bcftools query \
        --regions-file $snp_list_tsv \
        --samples-file $samples_file \
        --format "%CHROM\t%POS\t[%GT\t]\n" \
        --force-samples \
        $input_vcf \
        > $genotypes_GT_format


echo "Converting VCF genotypes to integers."
# First, write a row of column headers, giving "chr", "pos", then each line name
awk 'BEGIN{printf "chr\tpos"} {printf "\t%s", $0} END{print ""}' $samples_file > $genotypes_integer_format
# Convert GT format (0/0, 0/1, 1/0, 1/1) to integers (0,1,1,2)
# This is done for unphased and phased states (separated by / and | respectively.
cat ${genotypes_GT_format} | \
    sed \
        -e 's_0/0_0_g' \
        -e 's_0/1_1_g' \
        -e 's_1/0_1_g' \
        -e 's_1/1_2_g' \
        -e 's_0|0_0_g' \
        -e 's_0|1_1_g' \
        -e 's_1|0_1_g' \
        -e 's_1|1_2_g' \
        -e 's_./._NA_g' \
    >> $genotypes_integer_format


mkdir -p 05_results/16_gemma_ft12/output/06_top_SNP_genotypes
cp $genotypes_integer_format 05_results/16_gemma_ft12/output/06_top_SNP_genotypes