#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import pandas as pd
import h5py
import argparse
from warnings import warn
import itertools

in_options = argparse.ArgumentParser(
    description="""
    Compare offspring genotypes to an array of candidate parents.

    Divdes the genome into windows of roughly equal number of SNPs and estimates
    pi between the offsping and each candidate parent in each window. If the 
    offspring is sufficiently homozygous, pi should be close to zero for one of
    the parents in each window.

    To obtain a score for the validity of each cross, this identifies the minimum
    pi for either candidate at each window, take the log of one minus these 
    minima, and sums over each window. It returns a dataframe listing all possible
    pairs of combinations of parents ranked by score.    
    """)
in_options.add_argument(
    "--parental_geno", type=str, required=True,
    help="""
    Path to an HDF5 file giving genotype information on the panel of possible parents.
    This is usually created using `vcf_to_hdf5` from scikit-allel.
    """
)
in_options.add_argument(
    "--progeny_geno", type=str, required=True,
    help="""
    Path to an HDF5 file giving genotype information on the panel of possible parents.
    This is usually created using `vcf_to_hdf5` from scikit-allel.
    """
)
in_options.add_argument(
    "--progeny_name", type=str, required=True,
    help="Name of the progeny individual to test."
    )
in_options.add_argument(
    "--exp_parents", type=str, required=True, nargs=2,
    help="""
    The expected parents of the cross.
    Names should be strings separated by a space. For example:
        --exp_parents '6012' '9133'
    The parent IDs should be present in the genotype file for the parents.
    It is assumed that names in the genotype file are byte strings (e.g. b'6012'
    and b'9133'), but should be supplied here as regular strings.
    """
)
in_options.add_argument(
    "--snps_per_window", type=int, required=False, default=500,
    help="""
    Approximate number of SNPs in each window.
    The number of windows is the total number of SNPs on each chromosome divided
    by the number of SNPs per window, which is not always an integer.
    """
)
in_options.add_argument(
    "--outdir", type=str, required=True,
    help="Path to the directory to save the output files.")
args = in_options.parse_args()

# Load the HDF5 files for parents and progeny
parents_h5 = h5py.File(args.parental_geno, mode='r')
progeny_h5 = h5py.File(args.progeny_geno, mode='r')

# The index position of the progeny
progeny_ix = np.where(progeny_h5['samples'][:] == bytes(args.progeny_name, "utf8"))[0][0]
# Unique chromosome labels in the dataset
# Note: they are likely to be byte strings, and subsequent dictionary names will be as well.
chr_labels = np.unique(progeny_h5['variants']['CHROM'])

# array of names of the parents
parent_names = parents_h5['samples'][:]
# Check the names of the hypothesised parents are present in the genotype file.
for x in args.exp_parents:
    if bytes(x, 'utf-8') not in parent_names:
        raise ValueError(
            "One of the parents given ({}) is not present in the parental genotype file."
            )
        
# Dict to hold pi values for each chromosome
pi_genome_wide = {}

# Get pi across windows on each chromosome
for chr in chr_labels:
    # Get pi across windows on each chromosome
    # for chr in chr_labels:
    chr=chr_labels[0]
    # Integer index positions for the chromosome and SNP positions
    chr_ix = np.where(progeny_h5['variants']['CHROM'][:] == chr)[0]
    #Vector of indices defining where each window should begin
    start_points = np.arange(0, len(chr_ix), args.snps_per_window)
    # Empty dict to store the results for this chromosome
    pi_chr= {}

    for window in start_points[:5]:
        start = chr_ix[window]
        stop  = np.min([start + args.snps_per_window, chr_ix.max() ])

        # Arrays of genotypes in this window.
        # For the parents this has dimensions nloci x nparents x 2
        parents_window = parents_h5['calldata']['GT'][ start:stop ]
        # For the progeny this has dimensions nloci x 1 x 2
        progeny_window = progeny_h5['calldata']['GT'][ start:stop, progeny_ix]  
        progeny_window = np.expand_dims(progeny_window, axis=1)    
        # Use masked arrays to set negative values as invalid.
        parents_window = ma.masked_values(parents_window, -1)
        progeny_window = ma.masked_values(progeny_window, -1)
        
        # Calculate pi.
        array_of_differences = np.zeros(parents_window.shape[1])
        n_comparisons = np.zeros(parents_window.shape[1])
        # For each pair of parental and progeny haplotypes divide the number of
        # sequence differences by the number of loci that can be compared
        for i in [0,1]:
            for j in [0,1]:
                # Matrix of differences for parent haplotype i and progeny haplotype j
                diff_matrix = (parents_window[:,:,i] != progeny_window[:,:, j])              
                # Count the differences within the window, excluding NA sites
                array_of_differences += diff_matrix.sum(axis=0).data
                # Count the number of non-NA sites
                n_comparisons += ma.count(diff_matrix, axis=0)
                
        if any(array_of_differences > n_comparisons):
            raise ValueError("One or more comparisons in window starting at {} on chromosome {} has pi >= 1.".format(start, chr))
    
        # Calculate pi for this window.
        pi_chr[start] = array_of_differences / n_comparisons
    
    # Coerce to a dataframe with parental names as columns
    pi_chr = pd.DataFrame.from_dict(
        pi_chr, orient='index', columns=parent_names
    )
    pi_chr.insert(0, 'chr', chr)
    # send to genome0-wide dictionary
    pi_genome_wide[chr] = pi_chr
    
pi_genome_wide = pd.concat(pi_genome_wide)

# Score parentage
# Create a dataframe for each unique pair of candidates giving name and a score
# In each winow, similarity to either parent is 1-pi
# Total score is the sum of the (log) maximum of each pair.
parent_scores = []
for i,j in itertools.combinations(parent_names, 2):
    # Minimum pi in each window for either candidate
    max_values = (1-pi_genome_wide[[i,j]]).max(1)
    parent_scores = parent_scores + [{
        'parent1' : i.decode("utf-8") ,
        'parent2' : j.decode("utf-8") ,
        'score' : np.nansum(np.log(max_values))
    }]

# Sort by score
parent_scores = pd.DataFrame(parent_scores).sort_values('score', ascending=False)
parent_scores = parent_scores.iloc[range(10000)] # take only the top 1000 combinations

# Add a column indicating the candidates given as the true parents.
check_parent_AB = ((parent_scores['parent1'] == args.exp_parents[0]) & (parent_scores['parent2'] == args.exp_parents[1]))
check_parent_BA = ((parent_scores['parent2'] == args.exp_parents[0]) & (parent_scores['parent1'] == args.exp_parents[1]))
parent_scores['putative_parents'] = check_parent_AB | check_parent_BA
# write to disk.
parent_scores.to_csv(
    args.outdir + "/" + args.progeny_name + "_rank_candidates.csv",
    index=False
)

# Close the HDF5 files.
progeny_h5.close()
parents_h5.close()