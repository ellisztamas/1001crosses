# SNPs from the 1001 genomes project

File name:
`1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz`

Sequence data for the parents are taken from a VCF file produced by Fernando
Rabanal from

> Brachi, Benjamin, et al. "Plant genetic effects on microbial hubs impact host fitness in repeated field trials."Proceedings of the National Academy of Sciences 119.30 (2022): e2201285119

This is a reanalysis of previously published data from the [1001 genomes project](https://1001genomes.org/)
updated with newer SNP calling methods.

Note that accessions 1137 and 1074 are included as parents of crosses, but are
not in this VCF file.

# Imputed SNPs for the rest of the RegMap dataset

Filename:
`RegMap_panel/Arabidopsis_2029_Maf001_Filter80.vcf.gz`

From the publication:

> Arouisse, Bader, et al. "Imputation of 3 million SNPs in the Arabidopsis regional mapping population." The Plant Journal 102.4 (2020): 872-882.

Downloaded from here:
https://figshare.com/projects/Imputation_of_3_million_SNPs_in_the_Arabidopsis_regional_mapping_population/72887

# Resequenced data for 1435 and 5835

I suspected that these accessions are incorrect in the 1001 genomes data.
Seed stocks were resequenced, and data are in this folder, which is messy:

`/groups/nordborg/raw.data/athaliana/dna/the1001genomes_datasets/003.resequencing_seed_stocks/gbs_1001`

It is not easy to tell which samples are from which seed batch but I think I 
found one 1435 from 2017 and two from 5835 from the 2017 stock.
They are labelled

- `1435_seedstock.bam`
- `5835_seedstock_1.bam`
- `5835_seedstock_2.bam`

# Other files

- `segregating_parents.txt` is a list of parents I think are segregating somewhere in the genome.
- `admixture_groups_filiault_etal_2025.csv` gives the accessions and their 
    admixture groups for the lines used in:

> Brachi *et al*. (2025), Life-history trade-offs explain local adaptation in Arabidopsis thaliana, eLife14:RP107477