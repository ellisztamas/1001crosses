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
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=10
#SBATCH --time=08:00:00
#SBATCH --array=0-428

# Terminate the script automatically if any part of it fails (returns a non-zero exit status)
set -e

# Set working directory
source setup.sh

# Directory containing appropriately merged fastq files
indir=$workdir/03_merged_fastq
# Target directory for aligned reads
outdir=$workdir/04_aligned_bam
mkdir -p $outdir

# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas

# select files for read pairs 1 and 2.
read1_array=($(find $indir -type f -name '*val_1.fq.gz'))
read1_fastq=${read1_array[$SLURM_ARRAY_TASK_ID]}
# Swap read labels to get read2
read2_fastq=${read1_fastq/_R1_/_R2_}
read2_fastq=${read2_fastq/val_1/val_2}

# Name for the aligned bam file.
aligned_bam=${outdir}/$(basename -s _val_1.fq.gz $read1_fastq).sort.bam
aligned_bam=${aligned_bam/_R1_/_}
# Name for the deduplicated bam
dedup_bam=${aligned_bam/.sort/.sort.dedup}
dedup_metrics=${aligned_bam/.sort.bam/.sort.dedup.metrics.txt}
# Name for the BAM file with read-group information
read_group_bam=${aligned_bam/.sort/.sort.dedup.read_groups}




# Align reads to the genome
bwa mem -t 10 $genome $read1_fastq $read2_fastq | samtools sort -o $aligned_bam
if [ $? -eq 0 ] ; then echo "BWA completed successfully"; fi

# Mark duplicated reads
picard MarkDuplicates \
    I=$aligned_bam \
    O=$dedup_bam \
    M=$dedup_metrics
if [ $? -eq 0 ] ; then echo "Duplicate reads marked."; fi

# Add read-group information to the BAM file
# First, extract information on flow-cell, lane and library
flow_cell=$(samtools view $dedup_bam | head -n 1 | cut -d':' -f 3)
lane=$(samtools view $dedup_bam | head -n 1 | cut -d':' -f 4)
library=$(basename $dedup_bam | cut -d"_" -f 1)
# Add read groups to the BAM file
picard AddOrReplaceReadGroups \
    -I $dedup_bam \
    -O $read_group_bam \
    -RGID ${flow_cell}.${lane} \
    -RGLB ${library} \
    -RGPL ILLUMINA \
    -RGPU ${flow_cell}.${lane}.${library} \
    -RGSM $(basename $dedup_bam)
if [ $? -eq 0 ] ; then echo "Added read-group information successfully" ; fi

# Index bam files
samtools index $read_group_bam
if [ $? -eq 0 ] ; then echo "Indexing BAM files completed successfully"; fi