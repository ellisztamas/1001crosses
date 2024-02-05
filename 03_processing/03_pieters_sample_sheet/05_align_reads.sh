#!/usr/bin/env bash

# Align reads to the TAIR10 reference genome using BWA

# Input: trimmed, merged FASTQ files
# Output: BAM files aligned to the TAIR10 genome

# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023

# SLURM
#SBATCH --job-name=align_reads
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=10
#SBATCH --time=04:00:00
#SBATCH --array=0-428

# ml build-env/f2021
# ml bwa/0.7.17-gcc-10.2.0
# ml samtools/1.14-gcc-10.2.0

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
# i=$SLURM_ARRAY_TASK_ID
read1_array=($(find $indir -type f -name '*_R1_*.fastq.gz'))
read1_fastq=${read1_array[$SLURM_ARRAY_TASK_ID]}
read2_fastq=${read1_fastq/_R1_/_R2_}

# Name for the aligned bam file.
aligned_bam=${outdir}/$(basename -s .trim.fastq.gz $read1_fastq).sort.bam
aligned_bam=${aligned_bam/_R1_/_}
# Name for the deduplicated bam
dedup_bam=${aligned_bam/.sort/.sort.dedup}
dedup_metrics=${aligned_bam/.sort.bam/.sort.dedup.metrics.txt}



# Align reads to the genome
bwa mem -t 10 $genome $read1_fastq $read2_fastq | samtools sort -o $aligned_bam
if [ $? -eq 0 ] ; then echo "BWA completed successfully"; fi

# Mark duplicated reads
picard MarkDuplicates \
    I=$aligned_bam \
    O=$dedup_bam \
    M=$dedup_metrics
if[ $? -eq 0 ] ; then echo "Duplicate reads marked."; fi

# Index bam files
samtools index $dedup_bam
if [ $? -eq 0 ] ; then echo "Indexing BAM files completed successfully"; fi
