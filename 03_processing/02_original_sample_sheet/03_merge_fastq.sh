#!/usr/bin/env bash

# Merge FASTQ files that are split over lanes and run fastqc
#
# Samples from sequencing run H32TKDSX5 are split over two lanes.
# This script concatenates the split fastq files into a single files for each sample.
# It then runs fastqc and multiqc on the output
# 
# Input: 
    # Trimmed FASTQ file from two different sequencing runs, some split over two files
# Output:
    # A folder of FASTQ files merged where necessary
#
# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023

# SLURM
#SBATCH --job-name=merge_fastq
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=20GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00

# Set working directory
source setup.sh
# Location of unmerged fastq files
indir=${workdir}/02_trimmed_fastq
# Target directory for merged fastq files
outdir=${workdir}/03_merged_fastq
mkdir -p $outdir
# Directory for fastqc results
fastqc_dir=$outdir/qc
mkdir -p $fastqc_dir

files_to_merge=($(find $indir -type f -name '*_L003_*.trim.fastq.gz'))

# concatenate fastq files from run H32TKDSX5
# reads were split over two lanes and so need to be merged before alignement
for lane3 in ${files_to_merge[@]}
do
    lane4=${lane3/_L003_/_L004_}
    bs=$(basename $lane3)
    merged_file=${outdir}/${bs/_L003_/_}
    echo $merged_file
    cat $lane3 $lane4 > $merged_file
done

# copy files from HMN2MDRX2
# reads all went into the same lane, no need for merging
plate_5_files=($(find $indir -type f -name '219242*.trim.fastq.gz'))
for fastq in ${plate_5_files[@]}
do
    cp -v $fastq $outdir
done