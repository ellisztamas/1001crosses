#!/usr/bin/env bash
# 
# Unpack tarballs to a working directory.
# This uses the 'scratch-cbe' drive on the VBC CLIP cluster as the working 
# directory; change as necessary for your machine.

# Input: tarballs of raw fastq files
# OUtput: unpacked tarballs
#
# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023.

# SLURM
#SBATCH --job-name=unpack_F8s
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=20GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=0-2

i=${SLURM_ARRAY_TASK_ID}

# Path with tarballed data
indir=01_data/02_F8_unaligned_bams

# Set working directory
source setup.sh
# Output directory
outdir=$workdir/01_unzipped_fastq
mkdir -p $outdir

files=($indir/*tar.gz)

echo "File to unzip: ${files[$SLURM_ARRAY_TASK_ID]}"

FILE=${files[$SLURM_ARRAY_TASK_ID]}
if test -f "$FILE"; then
    echo "$FILE exists."
fi

# Unzip raw data to the working directory
tar -xzf ${files[$SLURM_ARRAY_TASK_ID]} --directory  ${outdir}