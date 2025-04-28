"""
Calculate the linkage disequilibrium (LD) between pairs of loci in an HDF5 file.

This script takes a HDF5 file containing genotype data and a text file
containing a list of SNPs that have been pruned for local LD. It calculates the
LD between pairs of loci in the HDF5 file that are in the target file, identifies
pairs with LD above a threshold and saves the results to a CSV file. SNPs are
processed in chunks to avoid memory issues.

Tom Ellis, 28th April 2025
"""
import numpy as np
import h5py
import allel
import pandas as pd
from time import time
import argparse

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
    '-r', '--threshold',
    help = "Threshold for LD. Pairs of loci with LD above this value will be saved.",
    type = float,
    default = 0.5
)
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
# snp_list=pd.read_csv('05_results/03_long_range_ld/output/parental_lines.prune.in', sep="\t", names=['snp'], dtype=str)
# args_threshold=0.3
# args_chunk_size=1000

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

print("Calculating LD between chunks of SNPs.")
# Empty dictionary to store the positions and LD values of pairs of SNPs
# that are above the linkage threshold
pairs_in_ld ={
    "i" : [], "j" : [], 'r2' : []
    }

t0=time()


i = args.chunk_ix
bin_i=bin_positions[i]
print(f"Processing chunk {i} of {len(bin_positions)}", flush=True)

# Subset the genotype calls from the HDF5 files for chunk i.
chunk_i = {}
chunk_i['calldata'] = geno_hdf5['calldata']['GT'][bin_i]
chunk_i['genotype_array'] = allel.GenotypeArray(chunk_i['calldata']).to_n_alt(fill=-1)

for j, bin_j in enumerate(bin_positions):
    if j >= i:
        # Subset the genotype calls from the HDF5 files for chunk j.
        chunk_j = {}
        chunk_j['calldata'] = geno_hdf5['calldata']['GT'][bin_j]
        chunk_j['genotype_array'] = allel.GenotypeArray(chunk_j['calldata']).to_n_alt(fill=-1)

        # Calculate the LD between the two chunks
        linkage_matrix_ij = allel.rogers_huff_r_between(chunk_i['genotype_array'], chunk_j['genotype_array'])**2
        linked_in_chunk_bool = linkage_matrix_ij > args.threshold

        # Only keep going if at least one pair of loci is above the threshold
        if linked_in_chunk_bool.any():
            # Indices of the pairs of loci that are above the linkage threshold
            # This is a two vectors giving the rows and columns of `calldata`
            linked_in_chunk_ix = np.where(linked_in_chunk_bool)

            # Use those indices to get the positions of the loci in the original
            # genotype data
            linked_in_hdf5_ix = {
                'i' : [ h5_chr_pos[locus] for locus in bin_positions[i][linked_in_chunk_ix[0]] ],
                'j' : [ h5_chr_pos[locus] for locus in bin_positions[j][linked_in_chunk_ix[1]] ]
            }

            # Store positions and LD values in `pairs_in_ld`
            pairs_in_ld = pd.DataFrame({
                'i' : linked_in_hdf5_ix['i'],
                'j' : linked_in_hdf5_ix['j'],
                'r2' : linkage_matrix_ij[linked_in_chunk_ix]
            })
            # Remove entries for LD with itself
            pairs_in_ld = pairs_in_ld.loc[pairs_in_ld['i'] != pairs_in_ld['j']]
            
            pairs_in_ld.to_csv(
                args.outfile,
                mode = 'a',
                header = False,
                float_format = '%.3f',
                index = False
            )
        
        else:
            print(f"Chunks {i} and {j} of {len(bin_positions)} has no pairs in LD", flush=True)
            continue

geno_hdf5.close()

print(f"Chunk {i} of {len(bin_positions)} took {time()-t0:.2f} seconds", flush=True)
