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
# Tom Ellis, 9th April 2025

# SLURM
#SBATCH --job-name=03_beagle
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=5GB
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=1-411

# Set working directory and load the conda environment
source setup.sh


# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# Table of valid F8 genotypes and their parents in .ped format
beagle_sample_sheet=$workdir/07_hmm_genotyping/01_sample_sheet/beagle_sample_sheet.ped

# Name of the F8 genotype to process
progeny_name=$(awk -F'\t' -v row=$i 'NR==row {print $2}' $beagle_sample_sheet)
# Names of the parents
parent1_name=$(awk -F'\t' -v row=$i 'NR==row {print $3}' $beagle_sample_sheet)
parent2_name=$(awk -F'\t' -v row=$i 'NR==row {print $4}' $beagle_sample_sheet)
# Value indicating whether the genotype is from the resequenced or low-coverage datasets
offspring_panel_name=$(awk -F'\t' -v row=$i 'NR==row {print $1}' $beagle_sample_sheet )

# VCF file containing all the SNPs from the published SNP matrix
reference_panel=$workdir/07_hmm_genotyping/02_reference_panel/reference_panel.vcf.gz
# VCF files for the progeny
# There are two, and the sample sheet has a column indicating which to use
resequenced_vcf=$workdir/08_resequencing/08_merge_VCF/resequenced.vcf.gz
low_coverage_vcf=$workdir/03_validate_genotypes/08_correct_split_vcf/F8_filtered.vcf.gz

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

echo "Processing progeny ${progeny_name} with parents ${parent1_name} and ${parent2_name}"

if [[ "$offspring_panel_name" == "resequenced" ]]; then 
    echo "Taking offspring data from the resequenced data set."
    offspring_panel=$resequenced_vcf
elif [[ "$offspring_panel_name" == "lowcoverage" ]]; then
    echo "Taking offspring data from the low-coverage data set."
    offspring_panel=$low_coverage_vcf
fi



# Parents
echo "Creating VCF file for the parents"
bcftools view \
    -s $parent1_name,$parent2_name \
    --output-type z \
    --output-file $cross_parents_unphased \
    $reference_panel
tabix $cross_parents_unphased

echo "Phasing the parents"
java -jar $beagle \
    gt=$cross_parents_unphased \
    out=$cross_parents_phased 
tabix ${cross_parents_phased}.vcf.gz



# Progeny
echo "Creating a VCF file for a single progeny individual."
bcftools view \
    -s $progeny_name \
    --output-type z \
    --output-file $single_progeny_unphased \
    $offspring_panel
tabix $single_progeny_unphased

echo "Get the row of the pedigree file for imputation."
awk -F'\t' -v row=$i 'NR==row' $beagle_sample_sheet > $ped_file

echo "Imputing the progeny."
java -jar $beagle \
    gt=$single_progeny_unphased \
    ped=$ped_file \
    ref=${cross_parents_phased}.vcf.gz \
    impute=true \
    impute-its=10 \
    out=$single_progeny_phased 
tabix ${single_progeny_phased}.vcf.gz
