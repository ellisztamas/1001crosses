#!/usr/bin/env bash

# Calculate lambda statistic for each GWAS result file
# 
# Run an Rscript to calculate lambda on each result file.
# Each script outputs a row giving cohort ID (parents, cohorts 1 or 2, or 
# permutd), replicate ID, number of loci, heritability and lambda.
# This row is appended to a text file.
#
# Tom Ellis, 21st July 2025

# SLURM
#SBATCH --job-name=07_calculate_lambda
#SBATCH --output=05_results/22_simulated_datasets/slurm/%x-%a.out
#SBATCH --error=05_results/22_simulated_datasets/slurm/%x-%a.err
#SBATCH --mem=1GB
#SBATCH --qos=rapid
#SBATCH --time=60:00
#SBATCH --array=1-799

# Set working directory and load the conda environment
source setup.sh

# === Inputs === 

i=${SLURM_ARRAY_TASK_ID}

# Array of results files for the parents only
input_array=($scratchdir/22_simulated_datasets/05_run_gemma/no_K/parents*assoc.txt)

# Pull out matching results files for the four cohorts
parents=${input_array[$i]}
cohort1=${parents/parents/cohort1}
cohort2=${parents/parents/cohort2}
permutd=${parents/parents/permutd}

# Path to the script to calculate lambda
lambda=05_results/22_simulated_datasets/06_calculate_lambda.R

# === Outputs ===

outdir=scratchdir/22_simulated_datasets/07_calculate_lambda
mkdir -p $outdir

# Temporary file to store the output of a single job.
tmp_out=$outdir/lambda_${i}.txt

# Final file to write to
outfile=05_results/22_simulated_datasets/output/lambda_values_noK.tsv


# === Main === 


# Fie headers giving the locus names and a column for r2
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo "This is the first task in the array. Creating output files."
    echo -e "cohort\trep\tnloci\th2\tlambda" > $outfile
else
    echo "This is not the first task in the array. Waiting for the first task to create output files."
    sleep ${SLURM_ARRAY_TASK_ID}
fi

# Run the script on each results file
echo "Running the script to calculate lambda"
Rscript $lambda --in $parents --out $tmp_out
Rscript $lambda --in $cohort1 --out $tmp_out
Rscript $lambda --in $cohort2 --out $tmp_out
Rscript $lambda --in $permutd --out $tmp_out

# Append the results file with these data.
echo "Writing to disk"
# flock -x "${outfile}.lock" -c "cat ${tmp_out} >> ${outfile}"
cat ${tmp_out} >> ${outfile}