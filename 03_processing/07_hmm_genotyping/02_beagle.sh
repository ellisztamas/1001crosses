#!/usr/bin/env bash

# Phase and impute individual F8 lines using Beagle.
#
# First pull out F8 genotypes identified as corrected by ibdpainting
# (see 03_processing/03_validate_genotypes) and extract the line name and names
# of the parents.
# Create a VCF file for the two parents only and phase them using Beagle.
# Then create a VCF file for the F8 line only and impute the missing genotypes
# using Beagle.
#
# Tom Ellis, 18th March 2024

# SLURM
#SBATCH --job-name=02_beagle
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=3GB
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=1-328

# Set working directory and load the conda environment
source setup.sh

# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# Table of valid F8 genotypes and their parents in .ped format
sample_list=03_processing/07_hmm_genotyping/output/sample_sheet.ped

# Name of the F8 genotype to process
progeny_name=$(awk -F'\t' -v row=$i 'NR==row {print $2}' $sample_list)
# Names of the parents
parent1_name=$(awk -F'\t' -v row=$i 'NR==row {print $3}' $sample_list)
parent2_name=$(awk -F'\t' -v row=$i 'NR==row {print $4}' $sample_list)

# VCF file containing all the SNPs from the published SNP matrix
all_parents=01_data/03_parental_genotypes/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz
# VCF file for the progeny with corrected sample names
# This is created by 03_processing/03_validate_genotypes/08_correct_split_VCFs.sh
all_progeny=$workdir/03_validate_genotypes/08_correct_split_vcf/F8_filtered.vcf.gz

# Path to the binary file for Beagle
beagle=02_library/beagle.r1399.jar

# === Output ===

outdir=$workdir/07_hmm_genotyping/${SLURM_JOB_NAME}
mkdir -p $outdir

# VCF file with only a single F8 and its two parents
# trio_vcf=$outdir/${progeny_name}.vcf.gz
# Pedigree file created by 01_phase_parents.sh
ped_file=$outdir/${progeny_name}.ped

# VCF files for a pair of parents 
cross_parents_unphased=$outdir/${progeny_name}_parents_unphased.vcf.gz
cross_parents_phased=$outdir/${progeny_name}_parents_phased
# VCF files for a single progeny
single_progeny_unphased=$outdir/${progeny_name}_unphased.vcf.gz
single_progeny_phased=$outdir/${progeny_name}_phased

# === Main ===

echo "Get the row of the pedigree file."
awk -F'\t' -v row=$i 'NR==row' $sample_list > $ped_file



# Parents
echo "Creating VCF file for the parents"
bcftools view \
    -s $parent1_name,$parent2_name \
    --output-type z \
    --output-file $cross_parents_unphased \
    $all_parents
tabix $cross_parents_unphased

echo "Phasing the parents"
java -jar $beagle \
    gl=$cross_parents_unphased \
    out=$cross_parents_phased 
tabix ${cross_parents_phased}.vcf.gz




# Progeny
echo "Creating a VCF file for a single progeny individual."
bcftools view \
    -s $progeny_name \
    --output-type z \
    --output-file $single_progeny_unphased \
    $all_progeny
tabix $single_progeny_unphased

echo "Imputing the progeny."
java -jar $beagle \
    gl=$single_progeny_unphased \
    ped=$ped_file \
    ref=${cross_parents_phased}.vcf.gz \
    impute=true \
    out=$single_progeny_phased 
tabix ${single_progeny_phased}.vcf.gz