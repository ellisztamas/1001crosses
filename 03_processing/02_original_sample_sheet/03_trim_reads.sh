#!/usr/bin/env bash

# Trim adaptors adapters from raw fastq files.

# Input: raw fastq files
# Output: The same fastq files without adapter sequences
#
# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023

# SLURM
#SBATCH --mem=20GB
#SBATCH --job-name=trim_adaptors
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=medium
#SBATCH --time=8:00:00
#SBATCH --array=0-812

i=$SLURM_ARRAY_TASK_ID

# DATA #
# Set working directory and load conda environment
source setup.sh
# Where the unzipped bam files are
indir=$workdir/01_unzipped_fastq
# Output file for trimmed files
outdir=$workdir/02_trimmed_fastq
mkdir -p $outdir

# Prepare file names
# Identify read pairs for matching fastq files (R1 and R2; I1 and I2 contain the barcodes)
read_array=($indir/**/demultiplexed/**/*_R1_*.fastq.gz)
read_pair_1=${read_array[$i]}
read_pair_2=${read_pair_1/_R1_/_R2_}
# prepare locations of trimmed fastq files
trimmed_read_1=${outdir}/$(basename -s .fastq.gz $read_pair_1).trim.fastq.gz 
trimmed_read_2=${outdir}/$(basename -s .fastq.gz $read_pair_2).trim.fastq.gz

# CUTADAPT #
cutadapt \
    -l 100 \
    -o $trimmed_read_1 -p $trimmed_read_2 \
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -q 20 \
    --minimum-length 20 \
    --pair-filter=any \
    $read_pair_1 $read_pair_2
if [ $? -eq 0 ] ; then echo "cutadapt completed successfully"; fi
