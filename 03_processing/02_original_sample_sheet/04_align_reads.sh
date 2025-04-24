#!/usr/bin/env bash

# Align reads to the TAIR10 reference genome using BWA, then mark duplicate reads
# add read-group information and index the resulting BAM files.
# 
# See the GATK page on read groups for an explanation of what this is doing:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
# This is needed to run base calibration with GATK later.
#
# Input:
    # trimmed, merged FASTQ files
# Output:
    # BAM files aligned to the TAIR10 genome, with duplicate reads marked and
        # read-group information.
    # .bai index files for each bam

# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023

# SLURM
#SBATCH --job-name=align_reads
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=10
#SBATCH --time=08:00:00
#SBATCH --array=0-428

# Terminate the script automatically if any part of it fails (returns a non-zero exit status)
set -e

# Set working directory
source setup.sh

# === Input ===

# Directory containing appropriately merged fastq files
indir=$scratchdir/03_merged_fastq
# Target directory for aligned reads
outdir=$scratchdir/04_aligned_bam
mkdir -p $outdir

# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas

# select files for read pairs 1 and 2.
read1_array=($(find $indir -type f -name '*val_1.fq.gz'))
read1_fastq=${read1_array[$SLURM_ARRAY_TASK_ID]}


# === Output ===

# Swap read labels to get read2
read2_fastq=${read1_fastq/_R1_/_R2_}
read2_fastq=${read2_fastq/val_1/val_2}

# Name for the aligned bam file.
aligned_bam=${outdir}/$(basename -s _val_1.fq.gz $read1_fastq).bam
aligned_bam=${aligned_bam/_R1_/_}
dedup_bam=${aligned_bam/.bam/.dedup.bam} 


# === Main ===

echo "Aligning reads to the genome."
bwa mem -t 10 $genome $read1_fastq $read2_fastq |\
    samtools sort -@ 10 -o $aligned_bam

echo "Marking and removing duplicates."
samtools fixmate -m $aligned_bam | \
    samtools markdup -Sro $dedup_bam

echo "Indexing bam files."
samtools index $dedup_bam