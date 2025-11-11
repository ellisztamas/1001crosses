#!/usr/bin/env bash

# Create samples of individual loci, stratified by minor-allele frequency (MAF).
# 
# This first subsets the parental VCF file by MAF, and saves marker names of 100
# random loci. It then that file to pull genotypes out of the parental and
# F8 VCF files. F8s are split by cohort.
#
# It returns a text file with a row for each marker and a column for each sample.
#
# Tom Ellis, 4th November 2025

# SLURM
#SBATCH --job-name=01_sample_genotypes
#SBATCH --output=05_results/24_single_locus_sims/slurm/%x-%a.out
#SBATCH --error=05_results/24_single_locus_sims/slurm/%x-%a.err
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-3

# Set working directory and load the conda environment
source setup.sh


# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# VCF files
# Parents
parents_vcf=03_processing/09_impute_haplotypes/output/parental_lines.vcf.gz
progeny_vcf=03_processing/09_impute_haplotypes/output/F8_imputed.vcf.gz

# Number of SNPs to subsample
nsnps=100

# Define lower and upper bounds on minor-allele frequencies.
maf_array=(0.0 0.1 0.2 0.3 0.4)
lower_maf=${maf_array[$i]}
upper_maf=${maf_array[$(($i+1))]}
range_maf="${lower_maf}-${upper_maf}"



# === Output ===

outdir=$scratchdir/24_single_locus_sims/01_sample_genotypes
mkdir -p $outdir

# Text files containing sample names
parents_names=$outdir/parents_names.txt
cohort1_names=$outdir/cohort1_names.txt
cohort2_names=$outdir/cohort2_names.txt

# A text file giving SNP positions to ensure that all datasets use the same markers
tmp_snplist=$outdir/tmp_snplist_${range_maf}.txt

# Prefix for output files without any file extension
parents_genotype_matrix=$outdir/parents_genotype_matrix_${range_maf}.tsv
cohort1_genotype_matrix=$outdir/cohort1_genotype_matrix_${range_maf}.tsv
cohort2_genotype_matrix=$outdir/cohort2_genotype_matrix_${range_maf}.tsv

# Project directory to stage out the results
projdir=05_results/24_single_locus_sims/output/01_samples_genotypes
mkdir -p $projdir



# === Main === 

# Get the names of crosses in each cohort to subset the VCF file later.
bcftools query -l $parents_vcf                > $parents_names
bcftools query -l $progeny_vcf | grep "_rep1" > $cohort1_names
bcftools query -l $progeny_vcf | grep "_rep2" > $cohort2_names


echo "Taking a random subsample of SNPs."
# bcftools view:
    # -q and -Q define a range of allele frequencies to retain.
    --min-ac: ensure there are at least two minor alleles
    -i 'F_MISSING<0.1': Ignore loci with <90% calls
# bcftools query exports a list of *all* SNPs matching the criteria above
# shuf takes a random sample of $nsnps of those snps.
bcftools view \
        -q "${lower_maf}:minor" -Q "${upper_maf}:minor" \
        --min-ac 2 \
        -i 'F_MISSING<0.1' \
        $parents_vcf | \
    bcftools query \
        -f "%CHROM\t%POS" | \
    shuf -n $nsnps \
    > $tmp_snplist



echo "Extracting genotypes for the parents."
# Write the first line as sample names.
echo -e "snp\t$(cat $parents_names | paste -sd \\t -)" > $parents_genotype_matrix
# Get the genotype calls
# bcftools query selects only the SNPs in $tmp_snplist and returns genotype calls at those sites
# awk converts VCF format (0/0, 0/1, 1/1, ./.) to integers (0,1,2,NA).
bcftools query \
        --regions-file $tmp_snplist \
        --format '%CHROM\_%POS[\t%GT]\n' \
        $parents_vcf |  \
    awk '{
        gsub(/0\/0|0\|0/, "0");
        gsub(/0\/1|1\/0|0\|1|1\|0/, "1");
        gsub(/1\/1|1\|1/, "2");
        gsub(/.\/.|.\|./, "NA");
        print
        }' \
    >> $parents_genotype_matrix



echo "Extracting genotypes for F8 cohort 1."
# Write the first line as sample names.
echo -e "snp\t$(cat $cohort1_names | paste -sd \\t -)" > $cohort1_genotype_matrix
# Get the genotype calls
bcftools query \
        --samples-file $cohort1_names \
        --regions-file $tmp_snplist \
        --format '%CHROM\_%POS[\t%GT]\n' \
        $progeny_vcf |  \
    awk '{
        gsub(/0\/0|0\|0/, "0");
        gsub(/0\/1|1\/0|0\|1|1\|0/, "1");
        gsub(/1\/1|1\|1/, "2");
        gsub(/.\/.|.\|./, "NA");
        print
        }' \
    >> $cohort1_genotype_matrix



echo "Extracting genotypes for F8 cohort 2."
# Write the first line as sample names.
echo -e "snp\t$(cat $cohort2_names | paste -sd \\t -)" > $cohort2_genotype_matrix
# Get the genotype calls
bcftools query \
        --samples-file $cohort2_names \
        --regions-file $tmp_snplist \
        --format '%CHROM\_%POS[\t%GT]\n' \
        $progeny_vcf |  \
    awk '{
        gsub(/0\/0|0\|0/, "0");
        gsub(/0\/1|1\/0|0\|1|1\|0/, "1");
        gsub(/1\/1|1\|1/, "2");
        gsub(/.\/.|.\|./, "NA");
        print
        }' \
    >> $cohort2_genotype_matrix


# Stage out the genotype files
cp $outdir/*_genotype_matrix_*.tsv $projdir