#!/usr/bin/env bash

usage () {
    cat <<HELP_USAGE

    A Bash script to combine a VCF and phenotype file to run through GEMMA.

    This subsets a genotype file to by the accessions in the phenotype file,
    creates the .bim, .bed, .fam and relatedness matrix files in Plink, and 
    passes these to GEMMA.

    For more on the input files, see the GEMMA manual at 
    http://www.xzlab.org/software/GEMMAmanual.pdf

    Tom Ellis, December 2023

    Options:
        --help        Print this help message.
        --vcf         VCF file containing SNP calls for the phenotyped lines.
                      If one or more lines are not present in the VCF file this 
                      will return and error.
        --pheno       Phenotype file formatted for GEMMA.
                      This should be tab delimited, usually with a column of zeroes,
                      a column of sample names and one or more columns of phenotypes.
        --outdir      Path to the directory to store the output.
        --covariates   Optional covariate file formatted for GEMMA.
                      This should be tab delimited with a column of ones to tell
                      GEMMA to fit an intercept, followed by one or more columns
                      of covariates.
        --gemma_args  Additional arguments to be passed to GEMMA.
                      Optional, but expect an error message if you don't specify
                      anything.

                      Note that for Bash-black-magic reasons this to work this
                      needs to be passed as a Bash variable in double quotes;
                      see the example.
    
    Output:
        Tab-delimited file giving the results of GEMMA with the same name as the
        phenotype file. Intermediate .bim, .bed and .fam files and a relatedness
        matrix are also created but are then deleted. Log files from PLINK and 
        GEMMA are kept.
    
    Example:
        # Paths to input files
        vcf=path/to/snp_matrix.vcf.gz
        pheno_file=path/to/phenotype_file.tsv
        covariates=path/to/covariate_file.tsv
        outdir=path/to/output_directory

        # Additional arguments need to be defined as a Bash variable in double 
        # rather than single quotes.
        gemma_args="-lmm 4 -maf 0.05"

        # Run the script
        02_library/run_GEMMA.sh --vcf \$vcf --phenotypes \$pheno_file --outdir \$outdir --covariates \$covariates --gemma_args "\${gemma_args}"
HELP_USAGE
}

# By default set no additional options to GEMMA
gemma_args=''

# *** Make sure you have a new enough getopt to handle long options (see the man page)
getopt -T &>/dev/null
if [[ $? -ne 4 ]]; then echo "Getopt is too old!" >&2 ; exit 1 ; fi

# Parse command line options using getopt
long_options="help,vcf:,,phenotypes:,outdir:,covariates::,gemma_args:"
options=$(getopt -o h --long "$long_options" -n "$(basename "$0")" -- "$@")
eval set -- "$options"

# Process command line options
while true; do
  case $1 in
    -h | --help)
      usage
      exit 0
      ;;
    --vcf) vcf=$2 ; shift 2;;
    --phenotypes) pheno_file=$2        ; shift 2;;
    --outdir) outdir=$2 ; shift 2 ;;
    --covariates) covariates=$2 ; shift 2 ;;
    --gemma_args)
      gemma_args="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    *)
      usage
      exit 1
      ;;
  esac
done

# Terminate the script automatically if any part of it fails (returns a non-zero exit status)
set -e

# Paths for output
tmp_dir=$outdir/tmp
with_K=$outdir/with_K
no_K=$outdir/no_K
mkdir -p $tmp_dir
mkdir -p $with_K
mkdir -p $no_K

# Pull out the basename of the phenotype file, without directories or the suffix
# This will be used to name output files.
pheno_name=$(basename $pheno_file)
pheno_name=${pheno_name%.*}

# Define temporary files
subsetted_VCF=$tmp_dir/${pheno_name}.vcf
plink_output=$tmp_dir/$pheno_name
relatedness_matrix=${pheno_name}_K_matrix

# Create PLINK file
echo "Creating the .bed .fam and .bim files, and filtering for samples with phenotypes..."
echo ""
plink2 \
  --vcf $vcf \
  --pheno $pheno_file \
  --make-bed \
  --set-missing-var-ids '#_@' \
  --keep $pheno_file \
  --out $plink_output

# Create relatedness matrix
echo "Calculating the relatedness matrix..."
gemma -bfile $plink_output \
  -gk 1 \
  -outdir $tmp_dir \
  -o $relatedness_matrix

# Run GEMMA
echo " Running GEMMA using the Plink relatedness matrix..."
gemma \
    -bfile $plink_output \
    -k ${tmp_dir}/${pheno_name}_K_matrix.cXX.txt \
    -c $covariates \
    -outdir $with_K -o $pheno_name \
    -lmm 4 \
    $gemma_args

echo " Running GEMMA without a relatedness matrix..."
gemma \
    -bfile $plink_output \
    -outdir $no_K \
    -o $pheno_name \
    -c $covariates \
    -lm 4 \
    $gemma_args

echo "Removing temporary files."
rm ${plink_output}.bim
rm ${plink_output}.bed
rm ${plink_output}.fam
rm ${plink_output}.log
rm ${plink_output}_K_*txt
