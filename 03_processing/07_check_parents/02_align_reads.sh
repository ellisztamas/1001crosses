#!/usr/bin/env bash

# Align reads to the TAIR10 reference genome using BWA, then mark duplicate reads

# Tom Ellis, 21st May 2025

# SLURM
#SBATCH --job-name=02_align_reads
#SBATCH --output=03_processing/07_check_parents/slurm/%x-%a.out
#SBATCH --error=03_processing/07_check_parents/slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=10
#SBATCH --time=08:00:00
#SBATCH --array=0-2

# Terminate the script automatically if any part of it fails (returns a non-zero exit status)
set -e

# Set working directory
source setup.sh

# === Input ===

i=$SLURM_ARRAY_TASK_ID

# Directory containing appropriately merged fastq files
indir=$scratchdir/07_check_parents/01_convert_and_trim

# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas

# select files for read pairs 1 and 2.
read1_array=($(find $indir -type f -name '*val_1.fq'))
read1_fastq=${read1_array[$i]}


# === Output ===

# Target directory for aligned reads
outdir=$scratchdir/07_check_parents/02_align_reads
mkdir -p $outdir

# Swap read labels to get read2
read2_fastq=${read1_fastq/_R1_/_R2_}
read2_fastq=${read2_fastq/val_1/val_2}

# Name for the aligned bam file.
aligned_bam=${outdir}/$(basename -s _val_1.fq.gz $read1_fastq).bam
aligned_bam=${aligned_bam/_R1_/_}
dedup_bam=${aligned_bam/.bam/.dedup.bam} 


# === Main ===

# echo "Aligning reads to the genome."
bwa mem -t 10 $genome $read1_fastq $read2_fastq > $aligned_bam

echo "Marking and removing duplicates."
samtools fixmate -@ 10 -m $aligned_bam - | \
    samtools sort -@ 10 | \
    samtools markdup -@ 10 -r - $dedup_bam

echo "Indexing bam files."
samtools index $dedup_bam