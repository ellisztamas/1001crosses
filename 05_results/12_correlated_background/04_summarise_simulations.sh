#!/usr/bin/env bash

# Run GEMMA on simulated datasets with and without the correction for the 
# relatedness matrix.
#
# Inputs:
#     VCF file for the parents
#     VCF file for the F8s
#     Text files with simulated phenotypes from 04_simulate_phenotypes.py
# Output:
#     Output files from GEMMA for each replicate

# Tom Ellis, 16th July 2024

# SLURM
#SBATCH --job-name=sum_corr_background
#SBATCH --output=05_results/12_correlated_background/slurm/%x-%a.out
#SBATCH --error=05_results/12_correlated_background/slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --array=0-1199

date 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source setup.sh

# === Paths to input files ===

i=$SLURM_ARRAY_TASK_ID

file_array=( $(cat /scratch-cbe/users/thomas.ellis/crosses/05_results/12_correlated_background/pheno_file_list.txt) )

parent_dir=$(dirname ${file_array[$i]})
rep1_dir=${parent_dir/parents/rep1}
rep2_dir=${parent_dir/parents/rep2}

echo $parent_dir

# === Output files === #

# Defined below

# === Script === #

# Summarise results
# Parents
python 02_library/summarise_sims.py -i $parent_dir/with_K/phenotype_file.assoc.txt > $parent_dir/with_K/sim_summary.txt
python 02_library/summarise_sims.py -i $parent_dir/no_K/phenotype_file.assoc.txt   > $parent_dir/no_K/sim_summary.txt

# Replicate 1
python 02_library/summarise_sims.py -i $rep1_dir/with_K/phenotype_file.assoc.txt > $rep1_dir/with_K/sim_summary.txt
python 02_library/summarise_sims.py -i $rep1_dir/no_K/phenotype_file.assoc.txt   > $rep1_dir/no_K/sim_summary.txt

# Replicate 2
python 02_library/summarise_sims.py -i $rep2_dir/with_K/phenotype_file.assoc.txt > $rep2_dir/with_K/sim_summary.txt
python 02_library/summarise_sims.py -i $rep2_dir/no_K/phenotype_file.assoc.txt   > $rep2_dir/no_K/sim_summary.txt





# Plot results
python 02_library/plot_gwas.py -i $parent_dir/with_K/phenotype_file.assoc.txt -o $parent_dir/with_K
python 02_library/plot_gwas.py -i $parent_dir/no_K/phenotype_file.assoc.txt   -o $parent_dir/no_K

# Replicate 1
python 02_library/plot_gwas.py -i $rep1_dir/with_K/phenotype_file.assoc.txt -o $rep1_dir/with_K
python 02_library/plot_gwas.py -i $rep1_dir/no_K/phenotype_file.assoc.txt   -o $rep1_dir/no_K

# Replicate 2
python 02_library/plot_gwas.py -i $rep2_dir/with_K/phenotype_file.assoc.txt -o $rep2_dir/with_K
python 02_library/plot_gwas.py -i $rep2_dir/no_K/phenotype_file.assoc.txt   -o $rep2_dir/no_K
