This directory contains scripts to infer recombination breakpoints in the F8
genotype data, assign ancestry to the haplotypes between them, and update 
F8 genotype data with parental data from the 1001-genomes project.

The idea is that:

1. Parents have good-quality genotype data
2. F8s have low coverage data

so if we know where regions of parental ancestry start and end we can just copy
that to the F8s.
Identifying the recombination breakpoints is done using a hidden-Markov model (HMM).
Rather than endlessly tweak the HMM to get it perfect, I instead wrote something
approximate and carefully checked each sample against the ibdpainting results.

Here is an overview of what the scripts do.
See the headers of each script for more detail.

* **01_sample_sheet.sh**: Create a sample sheet listing offspring, and which sequencing dataset to use for each.
* **02_identify_reliable_SNPs.sh**: Get positions of probably-syntenic SNPs
* **03_reference_panel.sh**: Create a VCF file for the parents at reliable SNP positions only.
* **04_emission_probabilities.R**: Function to get emission probabilities at each marker.
* **05_viterbi_hmm.R**: Function to get most likely ancestry at each marker.
* **06_run_HMM.R**: CLI script to run the HMM.
* **07_infer_breakpoints.sh**: SLURM submission script for the HMM.
* **08_check_breakpoints.\***: A report summarising changes required to breakpoint positions.
* **09_fill_haplotype_gaps.R**: Fill in the gaps in .bed files between parental haplotypes with F8 data.
* **10_splice_VCFs.sh**: Extract the regions of parental and offspring data for a single F8.
* **11_merge_VCF.sh**: Merge genotype data across samples.

Here are a list of samples that had poor or no data from ibdpainting results:

* 1435x6240_rep2
* 6019x6012_rep1
* 6024x1367_rep1
* 6030x992_rep1
* 6038x6039_rep2
* 6039x5835_rep1
* 6097x1318_rep2
* 6107x6128_rep1
* 6140x9388_rep1
* 6255x6136_rep1

In addition I also identified the following lines whilst compiling 
`08_check_breakpoints.qmd` that would benefit from better data:

* 6025x6046_rep2
* 5856x9452_rep2
* 6097x1318_rep2
* 6098x9451_rep1
* 6132x6255_rep2
* 6140x9388_rep1
* 6191x6965_rep2
* 9454x1363_rep1
* 9481x9399_rep1
* 9481x9399_rep1