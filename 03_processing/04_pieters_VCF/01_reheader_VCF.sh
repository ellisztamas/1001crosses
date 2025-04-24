#!/usr/bin/env bash

# Reheader the VCF for the F8s with sample names and split by replicate.

# Pieter's VCF file uses absolute paths to raw BAM files in the header
# This script pulls them out in order, replaces them with sample names, and 
# creates a new VCF file these shorter headers. Based on these names it then 
# creates separate VCF files for F8 replicates 1 and 2 (deprecated - this is
# done after imputing genotypes).

# Inputs:
#     Pieter's VCF file for the F8 plants
#     Pieter's calls on which samples could be validated (indirectly: this is done
#         in an R script, which this script calls)
# Outputs:
#     Zipped VCF files for the two replicates, excluding samples that could not
#         be validated, and with sample names as headers.
#     A bunch of intermediate VCF files that I should clean up
#     A tab-separated file giving variable positions among offspring.


# Tom Ellis, 15th December 2023

# SLURM
#SBATCH --job-name=reheader_VCF
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=10GB
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00

date 

# Load working directory and conda environment
source ./setup.sh
# Load R separately, since tidyverse doesn't play nicely with conda
ml build-env/f2022
ml r/4.2.0-foss-2021b

# == Inputs ==

# Path to the VCF file for the F8s created by Pieter
input_vcf=../crosses_continued/004.F8/001.genotyping/003.results/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC_crossesF8.vcf.gz

# === Output files ===

outdir=$scratchdir/04_pieters_VCF/01_reheader_VCF
mkdir -p $outdir

# Filename for the old headers
old_header=$outdir/vcf_header_to_change.txt
# Filename to save the new headers. This is created inside the R script, so don't mess with it.
new_header=$outdir/new_vcf_header.txt
# File with samples to be kept
samples_to_keep=$outdir/samples_to_keep.txt
# Path to save the resulting VCF file.
vcf_with_all_samples=$outdir/F8_snp_matrix.vcf.gz
# Path to save a VCF file with dubious samples removed
purged_VCF=$outdir/F8_snp_matrix_purged.vcf.gz
# Path after changing contig labels
chr_names=$outdir/chr_names.txt
vcf_with_Chr=$outdir/F8_snp_matrix_purged_chr.vcf.gz
# List of variable SNP positions
variable_positions=$outdir/variable_positions.tsv.gz

# === Script ===

# Pull out the existing header
bcftools query -l $input_vcf > $old_header
# Replace absolute paths with sample names, and save as $new_header
# Conda does not play nicely with r-tidyverse, so you may need to run this separately.
Rscript 03_processing/04_pieters_VCF/02_reorder_header.R

# # Swap the headers and keep only the samples confirmed by SNPmatch
echo "Swapping headers and removing dubious genotypes."
bcftools reheader --samples $new_header $input_vcf > $vcf_with_all_samples
bcftools view --samples-file $samples_to_keep $vcf_with_all_samples > $purged_VCF

# Swap chromosome labels from 1,2,3,4,5 to Chr1, Chr2, etc
# Create a space-delimited file for renaming chromosome labels in the VCF file.
echo "Renaming chromosome annotations."
cat > $chr_names << EOL
1 Chr1
2 Chr2
3 Chr3
4 Chr4
5 Chr5
EOL
# Rename chromosomes
bcftools annotate --rename-chrs $chr_names $purged_VCF | \
    bgzip -c > $vcf_with_Chr
tabix $vcf_with_Chr

# Split up the VCF file into 2
# Bash witchcraft to generate a comma-separated list of samples for each replicate:
# 1. extract sample names,
# 2. select only samples with 'rep1' in the name
# 3. Convert from having one sample per row to a comma-separated list.
# 4. Remove the trailing comma.
echo "Splitting the VCF file into replicates 1 and 2."
rep1_samples=$(bcftools query -l $vcf_with_Chr | grep "rep1" | tr '\n' ',' | sed -e 's/,$//')
rep2_samples=$(bcftools query -l $vcf_with_Chr | grep "rep2" | tr '\n' ',' | sed -e 's/,$//')
# Create separate VCF files for replicates 1 and 2
bcftools view -s $rep1_samples -Oz -o 03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.vcf.gz $vcf_with_Chr
bcftools view -s $rep2_samples -Oz -o 03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.vcf.gz $vcf_with_Chr
tabix 03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.vcf.gz
tabix 03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.vcf.gz

# Create a tab-separated file giving variable positions.
bcftools query -f"%CHROM\t%POS\n" $vcf_with_Chr | bgzip -c > $variable_positions

# # Tidy  up
# rm $old_header $new_header $samples_to_keep 
# rm $vcf_with_all_samples $purged_VCF
