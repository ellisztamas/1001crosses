"""
Create GEMMA-formatted phenotypes from SNP loci.

Pull out individual SNP loci from HDF5 files for the parents, and F8 data, and
turn these into phenotype files tha plink can read.

Note that this outputs to a working directory on scratch-cbe, which will only
work on he CLIP cluster! Normally scripts in this repo call setup.sh which 
defines a working directrory, but this won't work in Python out of the box.

Inputs:
    HDF5 files converted from VCF files for the parents and F8s.
    TSV file giving SNPs positions

Outputs:
    A TSV file for each SNP giving the genotype of each accession as a binary
        phenotype. Note: this is actually mean-centred and standardised to one
        standard deviation because plink converts zeroes to missing data.
        
Tom Ellis, 10th July 2024

"""

import h5py
import numpy as np
import pandas as pd
import os

# === Input files === #

# HDF5 versions of VCF files for the parents and F8s.
hdf5_files = {
    'parents' : "03_processing/04_pieters_VCF/output/parental_snp_matrix.hdf5",
    'rep1'    : "03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.hdf5",
    'rep2'    : "03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep2.hdf5"
}
geno = { k:h5py.File(path, mode = 'r') for k,path in hdf5_files.items() }

# TSV file with SNP positions
snps_to_simulate = pd.read_csv(
    "05_results/09_gwas_single_loci/output/snps_to_simulate.tsv", delimiter = "\t"
)

# === Output files === #


# outdir='/scratch-cbe/users/thomas.ellis/crosses/05_results/09_gwas_single_loci/'
outdir = '05_results/09_gwas_single_loci/output/tmp'

# === Script === #

for cohort in geno.keys():
    # Decode sample names from bytes to strings, and store for later.
    sample_names = [ x.decode('utf-8') for x in geno[cohort]['samples'][:] ]

    # Dictionary with an item for each chromosome.
    # Each item is a boolean list indicating whether a SNP in the Hdf5 file is on 
    # that chromosome. 
    # This is done once now because its really time consuming.
    chr_labels = snps_to_simulate['chrom'].unique()
    chr_ix = { chr : geno[cohort]['variants']['CHROM'][:] == str.encode(chr) for chr in chr_labels } 

    for i in range(snps_to_simulate.shape[0]):
        # Chromosome and position for this SNP.
        chr=snps_to_simulate.iloc[i].loc['chrom']
        pos=snps_to_simulate.iloc[i].loc['pos']
        snp_name = chr + "_" + str(pos)

        # Boolean list indicating position labels 
        pos_ix =  geno[cohort]['variants']['POS'][:] == pos
        # Find the position of the SNP indexing both chromosome and position.
        snp_ix = np.where(chr_ix[chr] * pos_ix)[0][0]

        # Create the phenotype file
        pheno_table = pd.DataFrame({
            'intercept'    : 0,
            'sample_names' : sample_names,
            'phenotype'    : geno[cohort]['calldata']['GT'][snp_ix].sum(1)
        })
        pheno_table = pheno_table.loc[pheno_table['phenotype'] >=0] # Remove missing data
        # Mean-centre and standardise to 1 st. deviation
        pheno_table['phenotype'] = pheno_table['phenotype'] - np.mean(pheno_table['phenotype'])
        pheno_table['phenotype'] = pheno_table['phenotype'] / np.std(pheno_table['phenotype'])

        # Write to disk
        write_dir = outdir + "/" + snp_name + "/" + cohort
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
        pheno_table.to_csv(
            write_dir + "/phenotype_file.tsv",
            sep = "\t",
            header=False,
            index = False
            )

{ v.close() for v in geno.values() }

