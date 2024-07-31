#!/usr/bin/env bash

# Compare offspring to its putative parents across the genome using check_trio.py
#
# Input:
    # HDF5 files for the parents and F8s prepared with scikit-allel.
    # VCF file for the F8s to extract sample names.
# Output:
    # Directory of CSV files giving pi between offspring and parent 1, offspring
    #   and parent2 and parent 1 and parent 2.
#
# Tom Ellis 9th January 2024

# SLURM
#SBATCH --job-name=check_trios
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-428

date

source setup.sh

# Input files
# HDF5 files for the parents and progeny
parents=03_processing/05_quality_control/output/hdf5_files/filtered_parental_SNP_matrix_mac20.h5
progeny=03_processing/05_quality_control/output/hdf5_files/F8_snp_matrix_mac20.h5
# Progeny VCF file (to extact offspring names)
progeny_vcf=03_processing/03_pieters_sample_sheet/output/F8_snp_matrix_mac20.vcf.gz

# parents=03_processing/05_quality_control/output/parents_test_data.h5
# progeny=03_processing/05_quality_control/output/F8_test_data.h5
# # Progeny VCF file (to extact offspring names)
# progeny_vcf=03_processing/05_quality_control/output/F8_test_data.vcf.gz

# Output directory
outdir=03_processing/05_quality_control/output/data_to_check_trios
mkdir -p $outdir

# Progeny sample name
sample_names=( $(bcftools query -l $progeny_vcf) ) # # Bash array of sample names
cross=${sample_names[$SLURM_ARRAY_TASK_ID]} # e.g. "9404x9390_F8_rep1"
# Get the names of the candidate parents for this
parent1="${cross%%x*}"  # Extracts everything before the 'x'
parent2=$(echo "${cross#*x}" | grep -oE '^[0-9]+')   # Extracts everything after the 'x' and before the first _

echo "Sample name: ${cross}"
echo "Parent 1: ${parent1}"
echo "Parent 2: ${parent2}"

# Run the script
02_library/check_trio.py \
  --parental_geno $parents \
  --progeny_geno $progeny \
  --progeny_name $cross \
  --snps_per_window 500 \
  --exp_parents $parent1 $parent2 \
  --outdir $outdir