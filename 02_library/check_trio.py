#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import pandas as pd
import h5py
import argparse
from warnings import warn

in_options = argparse.ArgumentParser(
    description="""
    Calculate similarity between and offspring and two candidate parents across
    their genomes

    Divdes the genome into windows of roughly equal number of SNPs and estimates
    pi between the offsping and each candidate parent in each window, and between
    the two parents. It returns a dataframe with a row for each window, giving
    pi between each pair of individuals.
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
# Index position of the parents in the array.
parents_ix = [ np.where(parent_names == bytes(x, "utf-8"))[0][0] for x in args.exp_parents ]

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
    # Integer index positions for the chromosome and SNP positions
    chr_ix = np.where(progeny_h5['variants']['CHROM'][:] == chr)[0]
    #Vector of indices defining where each window should begin
    start_points = np.arange(0, len(chr_ix), args.snps_per_window)
    # Empty dict to store the results for this chromosome
    pi_chr= {}

    for window in start_points:
        start = chr_ix[window]
        stop  = np.min([start + args.snps_per_window, chr_ix.max() ])
        
        # Arrays of genotypes in this window for the progeny and two parents
        # Dimensions: 3 individuals x nloci x 2 homologous chromosomes
        genotype_array = {
            'progeny' : progeny_h5['calldata']['GT'][start:stop, progeny_ix],
            'parent1' : parents_h5['calldata']['GT'][start:stop, parents_ix[0] ],
            'parent2' : parents_h5['calldata']['GT'][start:stop, parents_ix[1] ]
        }

        # Mark positions with NA data
        genotype_array = { k : ma.masked_values(v, -1) for k,v in genotype_array.items()}

        # Calculate pi.
        array_of_differences = np.zeros(3)
        n_comparisons = np.zeros(3)            
        # For each pair of parental and progeny haplotypes divide the number of
        # sequence differences by the number of loci that can be compared
        for i in [0,1]:
            for j in [0,1]:
                diff_dict = {
                    'progeny_vs_p1' : genotype_array['progeny'][:, i] != genotype_array['parent1'][:, j],
                    'progeny_vs_p2' : genotype_array['progeny'][:, i] != genotype_array['parent2'][:, j],
                    'p1_vs_p2'      : genotype_array['parent1'][:, i] != genotype_array['parent2'][:, j]
                }
                array_of_differences += np.array([ v.sum()     for v in diff_dict.values() ])
                n_comparisons        += np.array([ ma.count(v) for v in diff_dict.values() ])

        if any(array_of_differences > n_comparisons):
            raise ValueError("One or more comparisons in window starting at {} on chromosome {} has pi >= 1.".format(start, chr))

        # Calculate pi for this window.
        pi_chr[window] = array_of_differences / n_comparisons
        if any(pi_chr[window] >= 1):
            warn("One or more comparisons in window starting at {} on chromosome {} has pi >= 1.".format(start, chr))

    # Coerce to a dataframe with parental names as columns
    pi_chr = pd.DataFrame.from_dict(
        pi_chr, orient='index', columns=list(diff_dict.keys())
    )
    # send to genome-wide dictionary
    pi_genome_wide[chr] = pi_chr
    
# Reformat as a dataframe
pi_genome_wide = pd.concat(pi_genome_wide)
pi_genome_wide = pi_genome_wide.reset_index(names = ['chr', 'start']) # columns for chr and position
pi_genome_wide['chr'] = pi_genome_wide['chr'].str.decode("utf-8") # Chr labels as literal strings

pi_genome_wide.to_csv(
    args.outdir + "/" + args.progeny_name + "_check_trio.csv",
    index=False
)

# Close the HDF5 files.
progeny_h5.close()
parents_h5.close()