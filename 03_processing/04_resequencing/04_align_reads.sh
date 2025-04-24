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
#SBATCH --job-name=04_align_reads
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=15GB
#SBATCH --cpus-per-task=10
#SBATCH --time=02:00:00
#SBATCH --array=1-96

# Terminate the script automatically if any part of it fails (returns a non-zero exit status)
set -e

# Set working directory
source setup.sh

# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# Directory containing appropriately merged fastq files
indir=$scratchdir/04_resequencing/03_trim_reads

# Location of the reference genome to map to.
genome=01_data/01_reference_genome/TAIR10_chr_all.fas

# select files for read pairs 1 and 2.
read1_array=($(find $indir -type f -name '*val_1.fq.gz'))
read1_fastq=${read1_array[$i]}

# Swap read labels to get read2
read2_fastq=${read1_fastq/_R1_/_R2_}
read2_fastq=${read2_fastq/val_1/val_2}

# === Output ===

# Target directory for aligned reads
outdir=$scratchdir/04_resequencing/${SLURM_JOB_NAME}
mkdir -p $outdir

# Name for the aligned bam file.
aligned_bam=${outdir}/$(basename $read1_fastq)
aligned_bam=${aligned_bam/_val_1.fq.gz/_aligned.bam}
sorted_bam=${aligned_bam/aligned/sorted}
dedup_bam=${sorted_bam/sorted/dedup} 


# === Main ===

# echo "Aligning reads to the genome."
# bwa mem -t ${SLURM_CPUS_PER_TASK} $genome $read1_fastq $read2_fastq |\
#     samtools sort -@ ${SLURM_CPUS_PER_TASK} |\
#     samtools fixmate -@ ${SLURM_CPUS_PER_TASK} -m | \
#     samtools markdup -@ ${SLURM_CPUS_PER_TASK} -Sro $dedup_bam

# echo "Indexing bam files."
# samtools index $dedup_bam

echo "Aligning reads to the genome."
bwa mem \
    -t ${SLURM_CPUS_PER_TASK} \
    -o ${aligned_bam} \
    $genome \
    $read1_fastq $read2_fastq

# samtools sort \
#     -@ ${SLURM_CPUS_PER_TASK} \
#     -o $sorted_bam \
#     $aligned_bam

echo "Marking and removing duplicates."
samtools collate     -@ ${SLURM_CPUS_PER_TASK} -O -u ${aligned_bam} | \
    samtools fixmate -@ ${SLURM_CPUS_PER_TASK} -m -u - - | \
    samtools sort    -@ ${SLURM_CPUS_PER_TASK} -u - | \
    samtools markdup -@ ${SLURM_CPUS_PER_TASK} -r - $dedup_bam
# samtools fixmate -m $sorted_bam | \
#     samtools markdup \
#         --threads ${SLURM_CPUS_PER_TASK} \
#         --write-index \
#         -Sr \
#         $dedup_bam

# echo "Indexing bam files."
samtools index $dedup_bam