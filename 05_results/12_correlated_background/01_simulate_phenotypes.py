"""
Create PLINK-formatted phenotypes from SNP loci.

Pull out individual SNP loci from HDF5 files for the parents, and F8 data, and
turn these into case-control phenotype files thar plink can read. Simulate 
'focal SNPs' with a large effect on the phenotype, as well as a large number of
small effect loci across the genome.

Note that this outputs to a working directory on scratch-cbe, which will only
work on he CLIP cluster! Normally scripts in this repo call setup.sh which 
defines a working directrory, but this won't work in Python out of the box.

Inputs:
    HDF5 files converted from VCF files for the parents and F8s.
    TSV file giving SNPs positions

Outputs:
    A TSV file for each SNP giving the genotype of each accession as a binary
        phenotype. This is done for five different values for the main effect.
        Note: this is actually mean-centred and standardised to one
        standard deviation because plink converts zeroes to missing data.
        
Tom Ellis, 16th July 2024

"""

import h5py
import numpy as np
import pandas as pd
import os
import numpy.ma as ma
from scipy.special import expit


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

# Eigenvectors of a PCA on population structure
pop_structure = pd.read_csv(
    "05_results/01_pca/output/parental_snp_matrix.eigenvec", delimiter = " "
    )

# === Output files === #

outdir='/scratch-cbe/users/thomas.ellis/crosses/05_results/12_correlated_background'
# outdir = '05_results/12_correlated_background/output/tmp'

# === Script === #

# Dictionary with an item for each chromosome.
# Each item is a boolean list indicating whether a SNP in the Hdf5 file is on 
# that chromosome. 
# This is done once now because its really time consuming.
chr_labels = snps_to_simulate['chrom'].unique()
chr_ix = { chr : geno['parents']['variants']['CHROM'][:] == str.encode(chr) for chr in chr_labels } 

# Decode sample names from bytes to strings, and store for later.
sample_names={}
for cohort in geno.keys():
    sample_names[cohort] = [ x.decode('utf-8') for x in geno[cohort]['samples'][:] ]

# Randomly select SNP positions to serve as the genetic background
total_n_snps = geno['parents']['variants']['POS'].shape[0]
n_background_snps = 10000
background_snp_ix = np.sort(
    np.random.choice(
        np.arange(0, total_n_snps-1), size=n_background_snps, replace=False
        )
)
# Genotypes at the background SNPs, with missing data masked.
background_snps = {}
for cohort in geno.keys():
    background_snps[cohort] = ma.masked_less(
        geno[cohort]['calldata']['GT'][background_snp_ix].sum(2) -1,
        value = -1)

# Modify the sign of the allele labels so that the alternative allele is always
# correlated with PC1 of population structure
positive_cor_with_PC1 = [ np.corrcoef( i, pop_structure['PC1'])[0,1] > 0 for i in background_snps['parents'] ]
for cohort in geno.keys():
    background_snps[cohort] = background_snps[cohort] * np.where(positive_cor_with_PC1, 1, -1)[:, np.newaxis]


for i in range(snps_to_simulate.shape[0]):

    # Chromosome and position for the focal SNP.
    chr=snps_to_simulate.iloc[i].loc['chrom']
    pos=snps_to_simulate.iloc[i].loc['pos']
    snp_name = chr + "_" + str(pos)

    # Set an effect per locus such that positive alleles are associated with increased liability
    effect_per_locus = 0.0005

    for cohort in geno.keys():
        # Total background  effects on the logit scale
        logit_background = (effect_per_locus * background_snps[cohort]).sum(0).data

        # Boolean list indicating position labels 
        pos_ix =  geno[cohort]['variants']['POS'][:] == pos
        # Find the position of the SNP indexing both chromosome and position.
        snp_ix = np.where(chr_ix[chr] * pos_ix)[0][0]
        # Genotypes at the focal SNP
        focal_snp = geno[cohort]['calldata']['GT'][snp_ix].sum(1) -1
        focal_snp = ma.masked_less(focal_snp, value = -1)

        # Simulate five different effect sizes for the main effect.
        for p in [0.5, 0.6, 0.7, 0.8]:
            # Effect at the focal SNP, on the logit scale
            logit_p = np.log(p / (1-p)) * focal_snp
            # Overall liability is the (inverse logit) of the sum of effects at background and focal SNPs
            liability = expit(
                logit_p + logit_background
            )
            # Draw a binary phenotype. Note that this includes masked values
            binary_phenotype = np.random.binomial(1, liability)
            # Create the phenotype file
            pheno_table = pd.DataFrame({
                'intercept'    : 0,
                'sample_names' : sample_names[cohort],
                'phenotype'    : binary_phenotype + 1 # +1 so Plink recognises this as case-control data
            })
            # Remove rows with missing data at this locus
            if ma.is_masked(focal_snp):
                pheno_table = pheno_table.loc[~focal_snp.mask]
                
            # Write to disk
            write_dir = outdir + "/" + snp_name + "/" + cohort + "/p" + str(p)
            if not os.path.exists(write_dir):
                os.makedirs(write_dir)
            pheno_table.to_csv(
                write_dir + "/phenotype_file.tsv",
                sep = "\t",
                header=False,
                index = False
                )

{ v.close() for v in geno.values() }