# Reheader the VCF for the F8s with sample names

# Pieter's VCF file uses absolute paths to raw BAM files in the header
# This script pulls them out in order, replaces them with sample names, and 
# creates a new VCF file these shorter headers.

module load build-env/f2022
module load r/4.2.0-foss-2021b

# Path to the VCF file for the F8s created by Pieter
input_vcf=../crosses_continued/004.F8/001.genotyping/003.results/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC_crossesF8.vcf.gz
# Filename for the old headers
old_header=03_processing/pieters_VCF/vcf_header_to_change.txt
# Filename to save the new headers. This is created inside the R script, so don't mess with it.
new_header=03_processing/pieters_VCF/new_vcf_header.txt
# Path to save the resulting VCF file.
output_vcf=03_processing/pieters_VCF/F8_snp_matrix.vcf.gz

# Pull out the existing header
bcftools query -l $input_vcf > $old_header
# Replace absolute paths with sample names, and save as $new_header
Rscript 03_processing/pieters_VCF/01_reorder_header.R
# Swap the headers.
bcftools reheader --samples $new_header $input_vcf > $output_vcf

rm $old_header $new_header