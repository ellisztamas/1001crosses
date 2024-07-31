import pandas as pd
from glob import glob
import re
from numpy import median
import scipy.stats as st
import os

snp_directories = glob('05_results/09_gwas_single_loci/output/tmp/**/')

pvalue_list = []
gif_list = []

for i in range(len(snp_directories)):
    print(i)
    snp_id = re.search('Chr[1-5]_[0-9]+', snp_directories[i]).group(0)
    chr = re.search('Chr[1-5]', snp_id).group(0)
    pos = re.search('[0-9]+$', snp_id).group(0)

    pvalues = {
        'index' : i,
        'chr' : chr,
        'pos' : pos
    }
    gif = {
        'index' : i,
        'chr' : chr,
        'pos' : pos
    }


    for cohort in ['parents', 'rep1', 'rep2']:
        filename = snp_directories[i] + cohort + '/with_K/phenotype_file.assoc.txt'


        if not os.path.exists(filename):
            pvalues[cohort] = 'No file'
            gif[cohort]     = 'No file'

        else:
            result = pd.read_csv(filename, delimiter='\t')

            p_target_SNP = result.loc[(result['chr'] == int(chr[3])) & (result['ps'] == int(pos))]['p_lrt']
            if p_target_SNP.shape[0] == 1:
                p_target_SNP = p_target_SNP.to_numpy()[0]
            elif p_target_SNP.shape[0] == 0:
                p_target_SNP = "SNP missing"

            chi2 = st.chi2(df=1)
            lamb = chi2.isf(median(result.loc[~result['p_lrt'].isna()]['p_lrt'])) / chi2.median()

            pvalues[cohort] = p_target_SNP
            gif[cohort]     = lamb

    pvalue_list.append(pvalues)
    gif_list.append(gif)

pval_table = pd.DataFrame({
    'chr'     : [x['chr'] for x in pvalue_list],
    'pos'     : [x['pos'] for x in pvalue_list],
    'parents' : [x['parents'] for x in pvalue_list],
    'rep1'    : [x['rep1'] for x in pvalue_list],
    'rep2'    : [x['rep2'] for x in pvalue_list]
})

gif_table = pd.DataFrame({
    'chr'     : [x['chr'] for x in gif_list],
    'pos'     : [x['pos'] for x in gif_list],
    'parents' : [x['parents'] for x in gif_list],
    'rep1'    : [x['rep1'] for x in gif_list],
    'rep2'    : [x['rep2'] for x in gif_list]
})

pval_table.to_csv("05_results/09_gwas_single_loci/output/pvals_top_SNPs.tsv", sep="\t", index=False)
gif_table.to_csv("05_results/09_gwas_single_loci/output/genomic_inflation_factor.tsv", sep="\t", index=False)