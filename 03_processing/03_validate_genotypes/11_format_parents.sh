#!/usr/bin/env bash
#
# Subset the parental VCF file so it can be compared with the progeny.
#
# Subset the VCF with 1163 parents to include only those accessions that are the
# parent of at least one F8 line, and to include only those SNPs present in the 
# F8 VCF.
# 
# Tom Ellis, 6th March 2025

# SLURM
#SBATCH --job-name=format_parents
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --mem=5GB
#SBATCH --qos=rapid
#SBATCH --time=1:00:00

# Set a working directory and load the conda environment
source setup.sh

# === Input ===

parents_in=01_data/03_parental_genotypes/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz
progeny_in=$scratchdir/03_validate_genotypes/08_correct_split_vcf/F8_filtered.vcf.gz


# === Output ===

outdir=$scratchdir/03_validate_genotypes/11_format_parents
mkdir -p $outdir

# Data to subset parents
parent_names=$outdir/parent_names.txt
targets_file=$outdir/snps_to_keep.txt
# Output VCF file, filtered by sites and samples
parents_out=$outdir/parents_filtered.vcf.gz



# === Main === 

# Use the F8 line names to generate a list of unique parent names.
# This removes 1137 and 1074!
echo "Getting parental names."
bcftools query -l  $progeny_in | \
    awk -F'x' '{
        left = $1;
        right = $2;
        sub(/_rep[0-9]+$/, "", right);  # Remove "_rep" and any following digits from the right part
        print left;
        print right;
        }' | \
    sort | uniq | \
    grep -v "1137" | grep -v "1074" \
    > $parent_names

# create targets file
echo "Creating targets file."
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $progeny_in | bgzip -c > $targets_file
tabix -s1 -b2 -e2 $targets_file


echo "Subsetting the VCF file." 
bcftools view \
    -S $parent_names \
    -T $targets_file \
    -Oz \
    -o $parents_out \
    $parents_in
tabix $parents_out