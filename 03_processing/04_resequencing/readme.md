Scripts to process 91 F8 samples and two parents (1137 and 1074) that were
resequenced at higher coverage.

These samples were identified as:
- having haplotypes that didn't match the parents
- Having amiguous parentage
- having failed to sequence at all
- being a missing parent (accession 1074 and 1137)
- having a missing parent

Raw data is in `01_data/08_resequenced_F8s`.

These scripts call SNPs for these samples and validate the genotypes.

Lines that had no or poor data before, and still have no or poor data:
6030x992_rep1


Missing haplotype, but poor new data:


New missing haplotypes:
6036x1318_rep1
6105x5829_rep2
9450x9356_rep1
9450x9356_rep2

Thought to be missing haplotype, but actually correct
5831x8237_rep1
6092x6125_rep1
6104x6034_rep1
6109x9343_rep1
6115x1585_rep1
6918x9323_rep1
9353x6109_rep2

Thought to be missing haplotype, but actually wrong
5829x9336_rep1

Thought to be missing haplotypes, but poor or no new data
6019x6012_rep1
6019x6012_rep2
6038x6039_rep2
6097x1318_rep2
6098x9451_rep1
6107x6128_rep1
6136x6064_rep1
6136x6064_rep2
6140x9388_rep1
6255x6136_rep1

Thought to be missing haplotype, possible parent swap
9381x6198_rep1
9381x6198_rep2
9391x9353_rep1 # I thought this was correct the first time round
9391x9353_rep2

Cross appears incorrect
5829x9336_rep1
5829x9336_rep2
6030x992_rep2
8283x9471_rep1
8283x9471_rep2
9380x1063_rep1
9399x6097_rep1
9399x6097_rep2
