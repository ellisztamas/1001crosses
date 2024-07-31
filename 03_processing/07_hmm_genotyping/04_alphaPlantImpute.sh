#!/usr/bin/env bash

# Impute F8 genotypes using information from the parents using alphaPlantImpute2.
# 
# Inputs:
#     SNP matrix containing 215 parental accessions, subsetted from Fernando's 1163.        
#     SNP matrices for F8 replicates 1 and 2.
#     Text file giving each F8 and the expected parents for each.
# Outputs:
#     Plink files for the input SNP matrices
#     Library of haplotypes from the parents in .ped format.
#     Imputed genotypes, somehow.
#
# Steps:
#     Get a plantID for a single F8 genotype
#     Subset the offspring VCF file to contain only this genotype
#     Subset the parental VCF file to contain only the parents of the F8
#     Convert those two subsetted VCF files to plink format
#     Generate a library of haplotypes in the parents using alphaPlantImpute
#     Impute F8 genotypes using the haplotype library (returns a .ped file)
#     Convert the imputed .ped file back to VCF.
#     ????
#     Profit
#
# Tom Ellis, 28th June 2024

# SLURM
#SBATCH --job-name=alphaPlantImpute
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=0-399

date 
set -e

source ./setup.sh

# == Inputs ==

i=$SLURM_ARRAY_TASK_ID

# VCF files for the parents and F8s.
parental_input_vcf=03_processing/04_pieters_VCF/output/parental_snp_matrix.vcf.gz
offspring_input_vcf=03_processing/04_pieters_VCF/output/F8_snp_matrix_purged.vcf.gz

# Pull out the ID for a single F8 plant
offspring_sample_array=( $(bcftools query -l $offspring_input_vcf) ) # array of sample names
test_plantID=${offspring_sample_array[$i]}
echo "Running alphaPlantImpute script on sample ${test_plantID}."

# === Output files ===

outdir=$workdir/07_hmm_genotyping/04_alphaPlantImpute/${test_plantID}
mkdir -p $outdir


# PED format genotype files
parental_genotypes=$outdir/parental_genotypes
offspring_genotypes=$outdir/offspring_genotypes
# Text files listing parents and offspring in the input VCF files.
parents_list=$outdir/parents_list.txt
offspring_list=$outdir/offspring_list.txt

# Haplotype library created by alphaPlantImpute
haplotype_library=$outdir/haplotype_library
# Text file giving the parents
founders=$outdir/founders.txt
# Imputed genotype .ped file and .dosage from alphaPlantImpute
imputed_ped=$outdir/imputed

# Imputed VCF files
imputed_vcf_chr_as_integers=${outdir}/${test_plantID}_chr_as_integers
imputed_vcf_chr_as_strings=${outdir}/${test_plantID}.vcf.gz


# === Script ===

# Convert VCF file to PED
# (See also: http://apol1.blogspot.com/2014/11/best-practice-for-converting-vcf-files.html)
echo "Subsetting VCF files with BCFtools"
# Create a comma-separated list of the parent IDs to pass to BCFtools
# For example, change "6255x6136_rep1" to "6255,6136"
parents_list=${test_plantID/_rep[1,2]/}
parents_list=${parents_list/x/,}
bcftools view -s $parents_list $parental_input_vcf  -Ou -o ${parental_genotypes}.vcf
bcftools view -s $test_plantID $offspring_input_vcf -Ou -o ${offspring_genotypes}.vcf

# Prepare data files for alphaPlantImpute
echo "Converting VCF files to PED"
plink --vcf ${parental_genotypes}.vcf  --allow-extra-chr --const-fid --recode --out $parental_genotypes
plink --vcf ${offspring_genotypes}.vcf --allow-extra-chr --const-fid --recode --out $offspring_genotypes

# Create a file indicating the parents
# This is tab-delimited, with a column for the offspring, then parents. For example:
# 6255x6136_rep1 6255  6136
echo -e "${test_plantID}\\t${parents_list/,/\\t}" > $founders

# Run the program
echo "Creating the haplotype library."
AlphaPlantImpute2 -createlib \
    -maxthreads 8 \
    -out $haplotype_library \
    -ped $parental_genotypes.ped
echo "Imputing genotypes."
AlphaPlantImpute2 -impute \
    -out $imputed_ped \
    -library ${haplotype_library}.ped \
    -ped ${offspring_genotypes}.ped \
    -founder $founders


echo "Converting the imputed genotypes back to VCF."
# plink needs a map file with the same name as the .ped file, so make a copy
ln -f ${offspring_genotypes}.map ${imputed_ped}.map 
# Convert to (zipped) VCF
plink2 \
    --pedmap $imputed_ped \
    --recode vcf bgz \
    --out $imputed_vcf_chr_as_integers



# bcftools +fixref $imputed_vcf_chr_as_strings -- -f 01_data/01_reference_genome/TAIR10_chr_all.fas



# bcftools norm --check-ref e -f 01_data/01_reference_genome/TAIR10_chr_all.fas -Ou -o /dev/null $imputed_vcf_chr_as_strings



# bcftools +fixref \
#     -Ob -o /scratch-cbe/users/thomas.ellis/crosses/07_hmm_genotyping/04_alphaPlantImpute/6172x6148_rep2/6172x6148_rep2_fix.vcf.gz \
#     -- -f 01_data/01_reference_genome/TAIR10_chr_all.fas -m top \
#     $imputed_vcf_chr_as_string



# Swap chromosome labels from 1,2,3,4,5 to Chr1, Chr2, etc
# Create a space-delimited file for renaming chromosome labels in the VCF file.
echo "Renaming chromosome annotations."
cat > $outdir/chr_names.txt << EOL
1 Chr1
2 Chr2
3 Chr3
4 Chr4
5 Chr5
EOL
# Rename chromosomes and compress
bcftools annotate \
    --rename-chrs ${outdir}/chr_names.txt \
    -Oz -o ${imputed_vcf_chr_as_strings} \
    $imputed_vcf_chr_as_integers.vcf.gz

echo "Creating the index file for the final VCF file."
tabix ${imputed_vcf_chr_as_strings}