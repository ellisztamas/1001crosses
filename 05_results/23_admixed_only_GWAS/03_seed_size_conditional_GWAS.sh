#!/usr/bin/env bash

# Repeat the GWAS for seed size on lines descended from Northern and Southern 
# lines, but including the strongest association as a covariate.
#
# Tom Ellis, 13th October 2025

# SLURM
#SBATCH --job-name=03_seed_size_conditional_GWAS
#SBATCH --output=05_results/23_admixed_only_GWAS/slurm/%x-%a.out
#SBATCH --error=05_results/23_admixed_only_GWAS/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-1

# Load conda environment
source setup.sh 

# === Input files === #

i=$SLURM_ARRAY_TASK_ID

# Phenotype file
pheno_file_array=(05_results/23_admixed_only_GWAS/output/01_subset_lines/seed_size*tsv)
phenotype_file=${pheno_file_array[$i]}

# Pull out the basename of the phenotype file, without directories or the suffix
# This will be used to name output files.
pheno_name=$(basename $phenotype_file)
pheno_name=${pheno_name%.*}

# VCF files for parents and F8s
input_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz


# === Output files === #

# Output directory for GEMMA results and temporary files.
outdir=05_results/23_admixed_only_GWAS/output/03_seed_size_conditional_GWAS
mkdir -p $outdir

# Paths for output
tmp_dir=$outdir/tmp
with_K=$outdir/with_K
no_K=$outdir/no_K
mkdir -p $tmp_dir
mkdir -p $with_K
mkdir -p $no_K

# Define temporary files
subsetted_VCF=$tmp_dir/${pheno_name}.vcf
plink_output=$tmp_dir/$pheno_name
relatedness_matrix=${pheno_name}_K_matrix

# Text file sample names only.
sample_names=$outdir/tmp/${pheno_name}_sample_names.txt
# Text file with a column of ones, and a column of genotypes 
covariate_file=$outdir/tmp/${pheno_name}_covariate.txt

# File name for the output of GEMMA
gemma_file=${pheno_name}.assoc.txt

# Additional arguments passed to GEMMA
# -maf 0.05: Run all three kinds of statistical test using a minor-allele-frequency of 0.05
gemma_args="-maf 0.05"



# === Main === #

echo "VCF file: ${input_vcf}"
echo "Phenotype file: ${phenotype_file}"
echo "Output directory: ${outdir}"

# Create a covariate file with a column of ones, and a column of genotypes at the top SNP
# Comma-separated list of sample names
sample_names=$(cut -f2 "$phenotype_file" | paste -sd, -)
# Text file giving genotype calls at the top SNP, in sample order
# The sed command gives genotype calls as 0, 1, or 2 rather than 0/0, 0/1 or 1/1
# The awk command adds a column of ones, and makes the file tab-delimited.
bcftools query \
    --regions Chr5:16207115 \
    --samples $sample_names \
    --format "[%GT\n]" \
    $input_vcf | 
  sed 's/0\/0/0/; s/0|0/0/; s/0\/1/1/; s/1\/0/1/; s/0|1/1/; s/1|0/1/; s/1\/1/2/; s/1|1/2/; s/\.\/.*/NA/' | 
  awk -v OFS=' ' '{print 1, $0}' \
  > $covariate_file 


# # Create PLINK file
echo "Creating the .bed .fam and .bim files, and filtering for samples with phenotypes..."
echo ""
plink2 \
  --vcf $input_vcf \
  --pheno $phenotype_file \
  --make-bed \
  --set-all-var-ids '@_#' \
  --max-alleles 2 \
  --no-input-missing-phenotype \
  --keep $phenotype_file \
  --out $plink_output

# Create relatedness matrix
echo "Calculating the relatedness matrix..."
gemma -bfile $plink_output \
  -gk 1 \
  -outdir $tmp_dir \
  -o $relatedness_matrix

# Run GEMMA
echo " Running GEMMA using the Plink relatedness matrix..."
gemma \
    -bfile $plink_output \
    -k ${tmp_dir}/${pheno_name}_K_matrix.cXX.txt \
    -c $covariate_file \
    -outdir $with_K -o $pheno_name \
    -lmm 4 \
    $gemma_args

echo " Running GEMMA without a relatedness matrix..."
gemma \
    -bfile $plink_output \
    -outdir $no_K \
    -c $covariate_file \
    -o $pheno_name \
    -lm 4 \
    $gemma_args



# Create some plots
# Manhattan and QQ plots as .png files, which are quick
02_library/plot_gwas.py \
  --input $outdir/with_K/$gemma_file \
  --outDir $outdir/with_K
02_library/plot_gwas.py \
  --input $outdir/no_K/$gemma_file \
  --outDir $outdir/no_K

# Manhattan plots as an interactive HTML, which are slow
02_library/interactive_manhattan_plot.R \
    --input $outdir/with_K/$gemma_file \
    --threshold 1.5
02_library/interactive_manhattan_plot.R \
    --input $outdir/no_K/$gemma_file \
    --threshold 2
