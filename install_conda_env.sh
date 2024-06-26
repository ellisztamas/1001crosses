#!/usr/bin/env bash

# SLURM script to install the conda environment in the background
#
# Tom Ellis, 3rd Janurary 2024

# SLURM
#SBATCH --job-name=install_conda_env
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=20GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00

# module load anaconda3/2023.03
conda env create -f environment.yml