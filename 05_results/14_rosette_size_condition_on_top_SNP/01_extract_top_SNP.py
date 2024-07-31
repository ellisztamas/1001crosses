import pandas as pd
import numpy as np
import h5py


rep1 = pd.read_csv(
    "05_results/08_gemma_rosette_size/output/with_K/rosette_size_blups_F9_rep1.assoc.txt",
    delimiter = "\t"
    )
# Top SNP is Chr3:3738546
rep1.iloc[np.argmin(rep1['p_lrt'])]


path = "03_processing/04_pieters_VCF/output/F8_snp_matrix_purged_rep1.hdf5"
geno = h5py.File(path, mode = 'r')


ix = (geno['variants']['CHROM'][:] == str.encode('Chr3')) & (geno['variants']['POS'][:] == 3738546)

snp_table = pd.DataFrame({
    'sampleID' : [ x.decode('utf-8') for x in geno['samples'][:] ],
    'geno'     : geno['calldata']['GT'][ix].sum(2).squeeze()
    })

pheno = pd.read_csv(
    "03_processing/06_process_phenotypes/output/rosette_size_blups_F9_rep1.tsv",
    delimiter = "\t",
    names=['int', 'sampleID', 'phenotype']
    )


out = pd.merge(pheno, snp_table, how='left', on = "sampleID").\
    loc[:, ['int', 'geno']]

out['int'] += 1
out.loc[out['geno'] < 0, 'geno'] = -9

out.to_csv(
    "05_results/014_rosette_size_condition_on_top_SNP/output/chr3_3738546.tsv",
    index = False, header=False, sep = "\t"
    )    