# Permute genotypes at each locus for the parents.
#
# This applies the python script 01_permute_VCF.py to the VCF file for the 
# parents to randomise genotypes at each marker. This is done 200 times.
#
# Tom Ellis, 21st July 2025

# SLURM
#SBATCH --job-name=02_permute_VCF
#SBATCH --output=03_processing/08_simulated_datasets/slurm/%x-%a.out
#SBATCH --error=03_processing/08_simulated_datasets/slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=1-200

# Set working directory and load the conda environment
source setup.sh


# === Input ===

i=${SLURM_ARRAY_TASK_ID}

# VCF file to be permuted
in_vcf='03_processing/05_imputation/output/parental_lines.vcf.gz'


# === Output ===

outdir=$scratchdir/08_simulated_datasets/${SLURM_JOB_NAME}
mkdir -p $outdir
# Output VCF file with genotypes at each locus in a random order.
out_vcf=$outdir/permutd_$(printf "%03d\n" $i)


# === Main ===

# Permute genotypes
python 03_processing/08_simulated_datasets/01_permute_VCF.py \
    --input $in_vcf \
    --output $out_vcf \
    --seed $i
    
# Index the output file
tabix $out_vcf