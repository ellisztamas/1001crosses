"""
Create a table that will contain the information I can select SNPs from.

I am looking for SNPs of particular allele frequencies that are, or are not
correlated with population structure.

Input:
    HDF5 file containing genotype information on the parents.
    Table giving principle components of population structure from 05_results/01_pca
    Table of allele frequencies from 05_results/05_allele_freqs/

Returns:
    A tab-separated text file with a row for each SNP listing:
        Chromosome
        Position
        Number of non-missing alleles (total is sample size times 2)
        Minor allele frequency
        Correlation between SNP genotype and PC1
        Maximum correlation between SNP genotype and any of PCs 1 to 10
"""

import os
import h5py
import pandas as pd
import numpy as np
import numpy.ma as ma
import scipy.stats


# === Input files === #

# Allele frequencies estimated by vcftools
allele_freqs = pd.read_csv(
    "05_results/05_allele_freqs/output/parental_snp_matrix.frq",
    delimiter="\t",
    usecols=[0,1,3,5],
    header=0, names=['chrom', 'pos', 'mac', 'maf']
    )
# Change alternative-allele frequency to minor-allele frequency
allele_freqs['maf'] = np.where(allele_freqs['maf'] < 0.5, allele_freqs['maf'], 1-allele_freqs['maf'])

# Principle components of population structure
pca=pd.read_csv(
    "05_results/01_pca/output/parental_snp_matrix.eigenvec",
    delimiter=" "
)
# Genotype information
parents_HDF5_file = h5py.File(
    "03_processing/04_pieters_VCF/output/parental_snp_matrix.hdf5", 'r'
    )

# === Output file === #

outdir = '05_results/09_gwas_single_loci/output'
outfile = outdir + '/variant_table.csv'

if not os.path.exists(outdir):
    os.makedirs(outdir)

# === Script === #

n_SNPs = parents_HDF5_file['calldata']['GT'].shape[0]

cor_PC1=[]
max_cor=[]
mean_cor=[]

ix=range(0, n_SNPs, 10)

for i in ix:

    # Vector of genotypes at a single SNP, with missing data masked away
    masked_geno = ma.masked_array(
        parents_HDF5_file['calldata']['GT'][i],
        parents_HDF5_file['calldata']['GT'][i] < 0
        ).sum(1)
    # Correlation between genotypes at a single SNP and the first 10 PCs
    cor_with_PCs = np.array(
        [ scipy.stats.spearmanr(pca[col], masked_geno)[1] for col in pca.keys()[2:12] ]
    )
    # Correlation with PC1
    cor_PC1.append(cor_with_PCs[0])
    # Maximum correlation with any of the first 10 PCs
    max_cor.append(cor_with_PCs.max())
    # Mean correlation with any of the first 10 PCs
    mean_cor.append(cor_with_PCs.mean())


# Append columns and write to disk
allele_freqs.iloc[ix].\
    assign(cor_PC1 = cor_PC1).\
    assign(max_cor = max_cor).\
    assign(mean_cor = mean_cor).\
    to_csv(outfile, index=False)

parents_HDF5_file.close()