#!/usr/bin/env bash

# Assign variables common to all data processing scripts.
# If you need to change something, change it here.

# Tom Ellis, 27th November 2023

echo "Loading working directory and Conda environment from 00_setup.sh"

# Data processing was done on a special drive for large jobs on the VBC CLIP cluster
# This won't exist on other machines, so change it here if you have downloaded
# the code somewhere else.
workdir=/scratch-cbe/users/$(whoami)/crosses
export $workdir
# workdir=03_processing/output # example alternative working directory
echo "Working directory: ${workdir}" 
mkdir $workdir -p

# Load conda environment
# The first three lines are also CLIP-specific
# If you haven't already, install the environment with `conda env create -f environment.yml`
module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate 1001crosses
