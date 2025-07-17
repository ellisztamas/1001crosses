Scripts to check the identity of three parents that I think were incorrect
in the 1001 SNP matrix or are only in the RegMap panel:

- 1435
- 5835

We suspect that seeds for these accessions have been mixed up somewhere because:

- Meta data for these lines says they should be from Northern Sweden
- They actually look like inbred siblings from Southern Sweden
- The F8 crosses do not match the genotypes in the 1001-genomes SNP matrix
    - They do match the genomes in the RegMap panel

We think those used in the 1001 genomes project originated from the Bergelson lab.
We dont know if our stocks bulked in 2014 and 2017 are correct.
Fortunately these were sequenced at low coverage as part of Rahul Pisupati's
SNPmatch paper, so I can check.

Raw data from the resequencing effort are in
`/groups/nordborg/raw.data/athaliana/dna/the1001genomes_datasets/003.resequencing_seed_stocks/gbs_1001
`
Within that folder, filer for 1435 are:

- 2014 stock is somewhere in `06_C71EWANXX_20150609B_demux_2_rd4_1` but I don't 
know what the indices mean.
- 2017 stock: `03_C71EWANXX_20150609B_demux_3_rd3_1/C71EWANXX_3#708_505.bam`

and for 5835:

- 20_rd1_10_and_5_C9M8GANXX_20160806B_demux_7/C9M8GANXX_7#709_506.bam
- 13_rd2_10_and_rd1_n5_C89VWANXX_20160414B_demux_6/C89VWANXX_6#723_510.bam

I symlinked these to `01_data/03_parental_genotypes/`