"""
Calculate the linkage disequilibrium (LD) between pairs of loci in an HDF5 file.

This script takes a HDF5 file containing genotype data and a text file
containing a list of SNPs that have been pruned for local LD. It calculates the
LD between pairs of loci in the HDF5 file that are in the target file and saves
the results to a CSV file. SNPs are processed in chunks to avoid memory issues.

Tom Ellis, 28th April 2025
"""
import itertools
import numpy as np
import numpy.ma as ma
import h5py
import allel
import pandas as pd
from time import time
import argparse
from filelock import FileLock

# Parameters
parser = argparse.ArgumentParser(
    description = 'Identify pairs of loci in an HDF5 file that are in linkage disequilibrium (LD).')
parser.add_argument(
    '-g', '--genotypes',
    help = "Path to the HDF5 file containing the genotype data created by scikit allel.",
    required = True)
parser.add_argument(
    '-t', '--targets',
    help = "Path to a text file containing the SNPs to be used as targets for LD calculation, with one SNP per line. e.g. 'Chr1:123456'.",
    type = str,
    required = True)
parser.add_argument(
    '-c', '--chunk_size',
    help = "Size of the chunks to be used for LD calculation.",
    type = int,
    default = 1000)
parser.add_argument(
    '-i', '--chunk_ix',
    help = "Index of the chunk to be processed. This is used to split the SNPs into chunks for processing.",
    type = int, 
    required= True)
parser.add_argument(
    '-o', '--outfile',
    help = 'Path and filename for the output file. This will be a CSV file with the columns "i", "j", and "r2".',
    required = True)
args = parser.parse_args()



# === Inputs ===

# HDF5 file containing the genotype data
geno_hdf5 = h5py.File(args.genotypes, mode = 'r')

# Text file giving SNPs that have been pruned for local LD 
snp_list=pd.read_csv(args.targets, sep="\t", names=['snp'], dtype=str)

# parents_hdf5='03_processing/05_imputation/output/parental_lines.hdf5'
# geno_hdf5 = h5py.File(parents_hdf5, mode = 'r')
# snp_list=pd.read_csv('05_results/03_long_range_ld/output/parental_lines.snplist', sep="\t", names=['snp'], dtype=str)
# args.threshold=0.3
# args.chunk_size=5
# args.chunk_ix=2

# === Main ===

print("Loading SNP positions from HDF5 file that match the target file")
# Chr and position of the SNPs in the HDF5 file as lists of regular strings
h5_pos = geno_hdf5['variants']['POS'][:]
h5_pos = [ str(x) for x in h5_pos]
h5_chr = geno_hdf5['variants']['CHROM'][:]
h5_chr = [ x.decode('utf-8') for x in h5_chr ] # HDF5 stores strings as bytes
# Combine the chromosome and position into a single string
h5_chr_pos = [ x + ":" + y for x,y in zip(h5_chr, h5_pos) ]

# Get the indices of the SNPs in the HDF5 file that are in the target file
snps_in_target_file_ix = np.where(
    pd.Series(h5_chr_pos).isin(snp_list['snp'].values)
)[0]
# Split the indices into chunks of size chunk_size
bin_positions = np.array_split(
    snps_in_target_file_ix,
    np.floor(len(snps_in_target_file_ix) / args.chunk_size)
    )

def calculate_LD(chunk_i, chunk_j):
    """
    Calculate LD statistics between loci in two chunks of an HDF5 file.

    For each pair of loci this calculates D, D prime and r2.
    This ignores any locus which is missing or heterozygous in either sample
    
    It returns a dictionary of three arrays.
    """
    geno_array_i = allel.GenotypeArray(chunk_i).to_n_alt(fill=-9)
    geno_array_j = allel.GenotypeArray(chunk_j).to_n_alt(fill=-9)

    # Mask elements that are missing or heterozygous
    mask_i = (geno_array_i == -9) | (geno_array_i == 1)
    mask_j = (geno_array_j == -9) | (geno_array_j == 1)
    geno_array_i = ma.masked_array(geno_array_i, mask_i)
    geno_array_j = ma.masked_array(geno_array_j, mask_j)

    # Loop over pairs of loci and get LD statistics
    # The loop isn't ideal, but I am much more confident that this is doing what
    # it should than a vectorised solution. It is actually quicker than parsing
    # the locus names in the HDF5 file.
    d      = np.zeros([geno_array_i.shape[0], geno_array_j.shape[0]])
    dprime = np.zeros([geno_array_i.shape[0], geno_array_j.shape[0]])
    r2     = np.zeros([geno_array_i.shape[0], geno_array_j.shape[0]])
    for x in range(geno_array_i.shape[0]):
        for y in range(geno_array_j.shape[0]):
            # Dictionary of two-locus genotype frequencies
            # with keys 'p_00', 'p_02', 'p_20', 'p_22'
            pdict = {}
            for a in [0,2]:
                for b in [0,2]:
                    k = f'p_{a}{b}'
                    pdict[k] = np.nanmean((geno_array_i[x] == a) & (geno_array_j[y] == b))
            # Marginal frequencies of alt alleles at loci i and j
            pi = pdict['p_20'] + pdict['p_22']
            pj = pdict['p_02'] + pdict['p_22']
            # Expected and observed frequencies of the BB haplotype
            pipj = pi * pj
            pij  = pdict['p_22']

            # Calculate vanilla D            
            d_hat = pij - pipj
            d[x,y] = d_hat
            
            # If any allele is fixed, Dmax and r2 don't make sense
            if pi in [0,1] or pj in [0,1]:
                dprime[x,y] = np.nan
                r2[x,y]     = np.nan
            else:
                # Calculate D`
                negative_d = d_hat<=0
                dmax_if_negative = np.minimum(pi*pj, (1-pi)*(1-pj))
                dmax_if_positive = np.minimum( pi*(1-pj), (1-pi)*pj)
                dmax   = (negative_d * dmax_if_negative) + ((1-negative_d) * dmax_if_positive)
                dprime[x,y] = d_hat / dmax
                # Get r2
                r2[x,y] = (d_hat**2) / (pi * (1-pi) * pj * (1-pj))
    
    return {
        'd'      : d,
        'dprime' : dprime,
        'r2'     : r2
    }
    


print("Calculating LD between chunks of SNPs.")
t0=time()


i = args.chunk_ix
bin_i=bin_positions[i]
print(f"Processing chunk {i} of {len(bin_positions)}", flush=True)

# Subset the genotype calls from the HDF5 files for chunk i.
chunk_i = geno_hdf5['calldata']['GT'][bin_i]

for j, bin_j in enumerate(bin_positions):
    if j >= i:
        # Subset the genotype calls from the HDF5 files for chunk j.
        chunk_j = geno_hdf5['calldata']['GT'][bin_j]

        # Calculate the LD between the two chunks
        linkage_matrices_ij = calculate_LD(chunk_i, chunk_j)
        
        # Create a dataframe with all combinations of locus names in the two chunks.
        names_loci_i  = [ h5_chr_pos[locus] for locus in bin_positions[i] ]
        names_loci_j  = [ h5_chr_pos[locus] for locus in bin_positions[j] ]
        names_loci_ij = list(itertools.product(names_loci_i, names_loci_j))
        linkage_df_ij = pd.DataFrame(names_loci_ij, columns=['locus_i','locus_j'])

        # Flatten matrices of LD statistics and add as columns
        linkage_df_ij['d']      = [x for y in linkage_matrices_ij['d']      for x in y]
        linkage_df_ij['dprime'] = [x for y in linkage_matrices_ij['dprime'] for x in y]
        linkage_df_ij['r2']     = [x for y in linkage_matrices_ij['r2']     for x in y]

        # Append the results file.
        # FileLock is used to prevent multiple processes from writing to the file at the same time.
        with FileLock("output.csv.lock"):
            linkage_df_ij.to_csv(
                args.outfile,
                mode = 'a',
                header = False,
                float_format = '%.3f',
                index = False
        )

geno_hdf5.close()

print(f"Chunk {i} of {len(bin_positions)} took {time()-t0:.2f} seconds", flush=True)

