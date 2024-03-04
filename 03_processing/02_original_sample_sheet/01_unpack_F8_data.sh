#!/usr/bin/env bash
# 
# Unpack tarballs of raw sequence data for the F8s to a working directory.
#
# Input:
    # tarballs of raw fastq files
# Output:
    # unpacked tarballs
#
# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023.

# SLURM
#SBATCH --job-name=unpack_F8s
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=20GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array 0-2

# Set working directory
source setup.sh

# Path with tarballed data
indir=01_data/02_F8_unaligned_bams
# Output directories
outdir=$workdir/01_unzipped_fastq
mkdir -p $outdir
# Local project directory to copy the multiqc report
projdir=03_processing/02_original_sample_sheet/output/multiqc
mkdir -p $projdir


files=($indir/*tar.gz)
echo "File to unzip: ${files[$SLURM_ARRAY_TASK_ID]}"

infile=${files[$SLURM_ARRAY_TASK_ID]}
if test -f "$infile"; then
    echo "$infile exists."
fi

# Unzip raw data to the working directory
tar -xzf $infile --directory  ${outdir}



qc_dir=$outdir/$(basename -s .tar.gz $infile)/qc
cp ${qc_dir}/*multiqc*html $projdir