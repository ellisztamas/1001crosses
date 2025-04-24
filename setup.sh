#!/usr/bin/env bash

# Assign variables common to all data processing scripts.
# If you need to change something, change it here.

# Tom Ellis, 27th November 2023

echo "Loading working directory and Conda environment from setup.sh"

# Path to the project directory
# i.e. the directory where this script is located
projdir="$(dirname "$(readlink -f "$0")")"

# Data processing was done on a special drive for large jobs on the VBC CLIP cluster
# This won't exist on other machines, so change it here if you have downloaded
# the code somewhere else.
scratchdir=/scratch-cbe/users/$(whoami)/crosses
# scratchdir=03_processing/scratchdir # example alternative working directory

echo "Working directory: ${scratchdir}"
mkdir $scratchdir -p
if [ $? -ne 0 ] ; then 
echo "The working directory could not be created at the path given in setup.sh."
echo "Check that the directory exists on your machine and is writable."
fi 

# Load conda environment
# If you haven't already, install the environment with `conda env create -f environment.yml`
source activate 1001crosses