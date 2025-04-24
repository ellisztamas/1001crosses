#!/usr/bin/env bash
# 
# Unzip zipped sequencing data.
#
# Tom Ellis 20th March 2025.

# SLURM
#SBATCH --job-name=01_unzip_fastq
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --time=4:00:00


# Set working directory
source setup.sh

# === Input ===
# Path with tarballed data
infile=01_data/08_resequenced_F8s/22YHLCLT3_4_R18523_20250303.tar.gz



# === Output ===

# Output directories
outdir=$scratchdir/04_resequencing/${SLURM_JOB_NAME}
mkdir -p $outdir

# Local project directory to copy the multiqc report
projdir=03_processing/04_resequencing/output
mkdir -p $projdir



# === Main ===

# Unzip raw data to the working directory
tar -xzf $infile --directory  ${outdir}

cp $outdir/22YHLCLT3_4_R18523_20250303/qc/22YHLCLT3_4_R18523_multiqc_report.html $projdir