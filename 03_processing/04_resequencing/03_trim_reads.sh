#!/usr/bin/env bash

# Trim adaptors adapters from raw fastq files and run fastqc on the output.

# SLURM
#SBATCH --mem=2GB
#SBATCH --job-name=03_trim_reads
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=21,32,36,39,46,53,55,56,61,68,69,71



# Set working directory and load conda environment
source setup.sh

# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# Sample sheet giving BioSample,LibraryName,Run,refGenome,refPath,BioProject,fastq1,fastq2
sample_sheet=03_processing/04_resequencing/output/sample_sheet.csv


# === Output ===

# Output directory
outdir=$scratchdir/04_resequencing/${SLURM_JOB_NAME}
mkdir -p $outdir

# FastQC output directory
fastqc_dir=03_processing/04_resequencing/output/multiqc/
mkdir -p $fastqc_dir

# === Main ===

sample_name=$(awk -F',' -v row=$i 'NR==row {print $1}'  $sample_sheet)
fastq1=$(awk -F',' -v row=$i 'NR==row {print $7}'  $sample_sheet)
fastq2=$(awk -F',' -v row=$i 'NR==row {print $8}'  $sample_sheet)

# Trim Galore!
trim_galore \
    --paired \
    --nextera \
    --quality 20 \
    --length 20 \
    --clip_R1 15 --clip_R2 15 \
    --basename $sample_name \
    --fastqc \
    --trim-n \
    --output_dir ${outdir} \
    $fastq1 $fastq2

# Check if this is the last array task
if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]; then
    echo "Running fastqc on the trimmed reads"
    multiqc \
        --filename "multiqc_report_post_trimming.html" \
        --outdir $fastqc_dir \
        $outdir 
fi