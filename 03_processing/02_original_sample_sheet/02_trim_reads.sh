#!/usr/bin/env bash

# Trim adaptors adapters from raw fastq files and run fastqc on the output.
# 
# Input:
    # raw fastq files
# Output:
    # The same fastq files trimmed for adapter sequences, length, and leading
    # bases, and fastqc results on the output
#
# Tom Ellis, adapting code by Pieter Clauw, 24th November 2023

# SLURM
#SBATCH --mem=20GB
#SBATCH --job-name=trim_reads
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=0-812

# DATA #
# Set working directory and load conda environment
source setup.sh
# Where the unzipped bam files are
indir=$workdir/01_unzipped_fastq
# Output file for trimmed files
outdir=$workdir/02_trimmed_fastq
mkdir -p $outdir
# Output directory for fastqc results on the trimmed files
qcdir=$outdir/qc
mkdir -p $qcdir
# Output directory to save the multiqc report
multiqc_dir=03_processing/02_original_sample_sheet/output/multiqc/

# Prepare file names
# Identify read pairs for matching fastq files (R1 and R2; I1 and I2 contain the barcodes)
read_array=($indir/**/demultiplexed/**/*_R1_*.fastq.gz)
read_pair_1=${read_array[$SLURM_ARRAY_TASK_ID]}
read_pair_2=${read_pair_1/_R1_/_R2_}

# Trim Galore!
trim_galore \
    --paired \
    --nextera \
    --quality 20 \
    --length 20 \
    --clip_R1 15 --clip_R2 15 \
    --fastqc \
    --trim-n \
    --fastqc_args "--outdir ${qcdir}" \
    --output_dir ${outdir} \
    $read_pair_1 $read_pair_2
if [ $? -eq 0 ] ; then echo "Trim Galore! completed successfully"; fi

# If this is the last job pair of files to process, then also run multiqc
if [ $SLURM_ARRAY_TASK_ID -eq 858 ]
then
    multiqc -outdir ${multiqc_dir}/multiqic_after_trimming.html ${outdir}
fi