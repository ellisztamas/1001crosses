#!/usr/bin/env bash

# Assign variables common to all data processing scripts.
# If you need to change something, change it here.

# Tom Ellis, 27th November 2023

echo "Loading working directory and Conda environment from setup.sh"

# Data processing was done on a special drive for large jobs on the VBC CLIP cluster
# This won't exist on other machines, so change it here if you have downloaded
# the code somewhere else.
workdir=/scratch-cbe/users/$(whoami)/crosses
# workdir=03_processing/workdir # example alternative working directory
echo "Working directory: ${workdir}"
mkdir $workdir -p
if [ $? -ne 0 ] ; then 
echo "The working directory could not be created at the path given in setup.sh."
echo "Check that the directory exists on your machine and is writable."
fi 

# Load conda environment
# The first three lines are also CLIP-specific
# If you haven't already, install the environment with `conda env create -f environment.yml`
source activate 1001crosses

# Set an option to terminate any script when it encounters an error
set -e