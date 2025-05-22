import pandas as pd
from glob import glob
import re
import numpy as np
import scipy.stats as st
from statsmodels.stats import multitest
import os
import argparse

parser = argparse.ArgumentParser(
    description="Summarise results from a GEMMA analysis of simulated data. Prints to st. out."
)
parser.add_argument(
    '-i', '--input', help="A results file from a GEMMA linear model."
)
args = parser.parse_args()

filename = args.input

# Extract simulation information from the filename.
snp_id     = re.search('Chr[1-5]_[0-9]+', filename).group(0)
chr        = re.search('Chr[1-5]', snp_id).group(0)
pos        = re.search('[0-9]+$', snp_id).group(0)
cohort     = re.search('(parents|rep1|rep2)', filename).group(0)
liability  = re.search('p0\.[1-9]', filename).group(0)
includes_K = re.search('(with_K|no_K)', filename).group(0)

file_identifiers = [chr, pos, cohort, liability, includes_K]

if not os.path.exists(filename):
    print(",".join(file_identifiers, 'NA', "NA", "NA", "NA"))

else:
    result = pd.read_csv(filename, delimiter='\t')

    # Retreive the p-value and effect size for the focal SNP.
    ix = (result['chr'] == int(chr[3])) & (result['ps'] == int(pos)) # Row index
    p_target_SNP = result.loc[ix]['p_lrt']
    if p_target_SNP.shape[0] == 1:
        p_target_SNP = p_target_SNP.to_numpy()[0]
        beta = result.loc[ix]['beta'].to_numpy()[0]
    else:
        p_target_SNP = "NA"
        beta = "NA"
    file_identifiers.append(str(p_target_SNP))
    file_identifiers.append(str(beta))

    # Check for evidence of inflation at loci on *different* chromosomes
    # Subset those chromosomes, and filter NAs.
    other_chr_only = result.loc[
        (result['chr'] != int(chr[3])) & (~result['p_lrt'].isna())
        ]
    # Get the genomic inflation factor
    chi2 = st.chi2(df=1)
    lamb = chi2.isf(np.median(other_chr_only['p_lrt'])) / chi2.median()
    file_identifiers.append(str(lamb))
    # Count how many non-causal SNPs are above the Bonferroni threshold.
    bonferroni_cutoff = multitest.multipletests(result['p_lrt'], alpha = 0.05, method = 'fdr_bh')[3]
    false_pos_SNPs = (other_chr_only['p_lrt'] < bonferroni_cutoff).sum()

    file_identifiers.append(str(false_pos_SNPs))

    print(",".join(file_identifiers))