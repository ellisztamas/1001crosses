#!/usr/bin/env python


# import limix_plot as lp
from limix import plot
import pandas as pd
from numpy import log10
from statsmodels.stats import multitest
import math
import argparse
from os.path import basename

# Parameters
parser = argparse.ArgumentParser(
    description = 'Plot a Manhattan plot and QQ plot for a results file from GEMMA')
parser.add_argument(
    '-i', '--input',
    help = "Path to a tab-delimited results file from GEMMA. This should contain columns 'chr', 'pos' and 'p_lrt', meaning that GEMMA should be run with `-lmm 2` or `lmm -4`.",
    required = True)
parser.add_argument(
    '-o', '--outDir',
    help = 'Specify the output directory. All results will be saved in this directory.',
    required = True)
args = parser.parse_args()

# Import results file from GEMMA
gwas = pd.read_csv(
    args.input, delimiter="\t",
    usecols = ['chr', 'ps', 'p_lrt']
    )
gwas = gwas.rename(columns={
    "chr" :"chrom",
    "ps": 'pos',
    'p_lrt': "pv"
    })

# Name of the file, without directory path or file extension.
file_prefix = basename(args.input).replace('.assoc.txt', '')

# Calculate the Bonferroni-adjusted significance threshold.
bonferroni_cutoff = 0.05 / gwas.shape[0]

Bonferroni = multitest.multipletests(gwas['pv'], alpha = 0.05, method = 'fdr_bh')[3]
# # Plot the Manhattan plot
# lp.manhattan(gwas)
# plt = lp.get_pyplot()
# _ = plt.axhline(-log10(bonferroni_cutoff), color='red')
# _ = plt.title(file_prefix)
# plt.savefig(args.outDir + "/" + file_prefix + "_manhattan_plot.png")
# plt.close()

# # Plot a QQ plot
# lp.qqplot(gwas['pv'])
# plt=lp.get_pyplot()
# plt.title(file_prefix)
# plt .savefig(args.outDir + "/" + file_prefix + "_qqplot.png")

# plot results
# Manhattan plot
plot.manhattan(gwas)
plt = plot.get_pyplot() 
_ = plt.axhline(-math.log(Bonferroni, 10), color='red')  
plt.savefig(args.outDir + "/" + file_prefix + "_manhattan_plot.png")
plt.close()

# QQ-plot
plot.qqplot(gwas['pv'])
plt = plot.get_pyplot()
plt.savefig(args.outDir + "/" + file_prefix + "_QQplot.png")
plt.close()