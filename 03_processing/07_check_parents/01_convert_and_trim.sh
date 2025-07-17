#!/usr/bin/env bash

# Split raw bam files by read-mate pair and convert from BAM to fastq.
# Trim adaptors adapters from fastq files and run fastqc on the output.
#
# Tom Ellis, 21st May 2025

# SLURM
#SBATCH --mem=10GB
#SBATCH --job-name=convert_and_trim
#SBATCH --output=03_processing/07_check_parents/slurm/%x-%a.out
#SBATCH --error=03_processing/07_check_parents/slurm/%x-%a.err
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=0-2

source setup.sh

# === Input ===

i=$SLURM_ARRAY_TASK_ID

# Directory with the BAM files to be genotyped
bam_array=(01_data/03_parental_genotypes/*bam)
raw_bam=${bam_array[$i]}


# === Output ===

# Directory to store the output files
outdir=$scratchdir/07_check_parents/01_convert_and_trim
mkdir -p $outdir

# Sorted BAM file
sorted_bam=$outdir/$(basename $raw_bam)
# Fastq files separated by read pair
fastq_1=${sorted_bam/.bam/_R1.fastq}
fastq_2=${sorted_bam/.bam/_R2.fastq}

# Output directory for fastqc results on the trimmed files
qcdir=$outdir/qc
mkdir -p $qcdir



# === Main ===

echo "Converting BAM to fastq."
samtools sort -n -o $sorted_bam $raw_bam
bedtools bamtofastq -i $sorted_bam -fq $fastq_1 -fq2 $fastq_2

echo "Running trim galore on the fastq file."
trim_galore \
    --paired \
    --quality 20 \
    --length 20 \
    --clip_R1 15 --clip_R2 15 \
    --fastqc \
    --trim-n \
    --fastqc_args "--outdir ${qcdir}" \
    --output_dir ${outdir} \
    $fastq_1 $fastq_2
