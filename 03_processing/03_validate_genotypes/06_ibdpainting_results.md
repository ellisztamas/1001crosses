
# Validation of the ancestry of F8 genotypes

This document gives an overview of visual valdidation of genotypes from the 
1001 crosses project using the program [ibdpainting](https://github.com/ellisztamas/ibdpainting).
It gives representative examples of different outcomes and a summary of how many cases were found of each.
A comprehensive list of the result for each line is given in the file `ibdpainting_results.csv`

## Contents

* [Introduction](#introduction)
    - [Crossing experiment](#crossing-experiment)
    - [IBD painting](#ibd-painting)
* [Samples that worked](#samples-that-worked)
    - [Matching both parents](#matching-both-parents)
    - [Residual heterozygosity](#residual-heterozygosity)
* [Cases that are easy to diagnose](#cases-that-are-easy-to-diagnose)
    - [No data](#no-data)
    - [Expected parents are not in the reference panel](#expected-parents-are-not-in-the-reference-panel)
    - [Labels were swapped](#labels-were-swapped)
    - [Offspring is actually self-pollinated](#offspring-is-actually-self-pollinated)
+ [Ambiguous cases](#ambiguous-cases)
    - [Missing haplotypes](#missing-haplotypes)
    - [Parent appears incorrect](#parent-appears-incorrect)
    - [One parent *might* have been swapped](#one-parent-might-have-been-swapped)
* [Summary and action](#summary-and-action)

## Background

#### Crossing experiment

Random pairs of 219 accessions of *Arabidopsis thaliana* from Sweden were crossed, and one plant from each cross was then maintained by self-pollination.
For most crosses two F2 seeds were propogated to give two replicates per cross.
At the F8 generation 429 lines were genotyped with Illumina short-read data at 
about 6X coverage.
**We want to confirm that the F8 genotypes really match the parents 
they ought to come from**.

To check this, we calculate genetic distance of each F8 to to the two expected
parents in windows of 500kb using the program [ibdpainting](#ibd-painting).
This is done for SNPs inside genes only, because SNP calls outside genes are dubious due to widespread structural polymorphism.
We also calculate the distance to each of 1163 accessions based on data from an extended set of the 1135 accessions from the 1001 genomes project, from this paper:

> Brachi, Benjamin, et al. "Plant genetic effects on microbial hubs impact host fitness in repeated field trials." Proceedings of the National Academy of Sciences 119.30 (2022): e2201285119

One goal would be to determine whether parents are correct so that we can use the high-quality parental genotypes to improve the low-coverage F8 genotyping data.

#### IBD painting

[ibdpainting](https://github.com/ellisztamas/ibdpainting) is a visual tool to validate ancestry across a genome.
For windows across the genome this calculates the genetic distance from a test individual to every genotype in a reference panel.
It then plots these distances for the expected parents, and a subsample of the most likely parents.
Smaller distances mean a closer match.

In addition, it also generates an IBD score file for each *pair* of candidate parents.
In each window it takes the minimum distances to either parent, and takes the average distance over the whole genome.

## Samples that worked

#### Matching both parents

F8 individuals have been through seven rounds of self-fertilisation, so we expect them to be homozygous across the genome, with megabase-sized haplotypes matching one of the two parents.
Figure 1 shows a good example of how such a genome should look.
Each point shows the average genetic distance from the sample 1_A3 in each window to the two expected parents (8423 and 9413 in red and blue).
Also shown in grey are distances to the ten next-best-matching accessions.
You can see that the F8 genome is a perfect match for one of the two expected parents at each window in the genome (where distance is zero), and that there are megabase-scale blocks of identity.
We can also see that distances are not quite zero where ancestry blocks swap over; that's fine, it just means that there is a recombination site somewhere in the 500kb window, such that neither parent is a perfect match.
A bit more than half of the samples look like this.

![Figure 1: An example of an F8 individual which is homozygous for one of the two parents at each window in the genome.](output/ibdpainting/1_A3_plot_ibd.png)

Figure 2 shows a similar example of an F8 that matches one or both parents in each window.
In this case however the distances are not zero as they were in figure 1, and are generally noisy.
The notes from the library preparation reveal that this sample had a very low DNA concentration.
Nevertheless, it does match the expected parents, so we can say that the seed material is likely not contaminated, although one might wish to repeat the genotyping.
There are seven samples that look like this.

![Figure 2: An example of an F8 individual which matches one parent across the genome, but the signal is noisy](output/ibdpainting/4_F9_plot_ibd.png)

#### Residual heterozygosity

Figure 3 shows an example of an F8 individual that matches one of the two expected parents, except for a two regions on chromosome 5.
Nevertheless, the expected parents are still the closest matches.
A simple explanation is that the F8 is heterozygous for two parental haplotypes in these regions.
This can be confirmed by checking the genotypes in the region with `bcftools`:
```
bcftools query -s 1_C12 -r Chr5:12000000-13000000 -f"%CHROM %POS [ %GT] [ %AD]\n" vcf_file.vcf.gz | grep "0/1" | wc -l
```
In this collection, genuinely heterozygous regions like this one have hundreds of heterozygous SNP calls per megabase.
Homozygous regions have tens of heterozygous SNP calls per megabase, which are almost certainly miscalled.

![Figure 3: An example of an F8 individual with residual heterozygosity on chromosome 5.](output/ibdpainting/1_C12_plot_ibd.png)

There are 53 samples that show clear signs of residual heterozygosity in at least one place in the genome.
This is a lower bound on the true number, because it is difficult to spot heterozygosity in cases where a parent is missing, or genotyping is ambiguous (see cases below).

## Cases that are easy to diagnose

#### No data

16 samples are recorded as having zero or very low concentrations of DNA in the
NGS libraries.
Not surprisingly, it was not possible to call SNPs at almost any window, so the plots look nearly empty (figure 4).
This should be easy to fix by sowing more seeds and repeating the genotyping.

![Figure 4: An example of an F8 individual with residual heterozygosity on chromosome 5.](output/ibdpainting/1_E7_plot_ibd.png)

#### Expected parents are not in the reference panel

Figure 5 shows a case where nothing matches.
Lines for the expected parents are not plotted at all.
The figure shows only grey lines, none of which match the F8.
It turns out that the expected parents (accessions 1137 and 1074) were not included in the reference panel.

![Figure 5: The parents of this cross were not in the reference panel, so no blue or red lines are plotted at all.](output/ibdpainting/3_D5_plot_ibd.png)

Likewise, figure 6 shows an F8 from a cross between accessions 1137 and 8241.
You can see long tracts where 8241 matches the F8, and long tracts where nothing matches, which are presumably homozygous for the 1137 haplotypes.

![Figure 6: An example where one parent is missing from the reference panel, but the other is present, and a good match.](output/ibdpainting/1_D8_plot_ibd.png)

Four crosses are affected by the two missing accessions (including the replicates of each, that gives eight lines in total)

* 1137x1074
* 8241x1137
* 1137x6237
* 6195x1074

#### Labels were swapped

Sample 2_H4 doesn't match it expected parents 9481 and 8307 (figure 7).
However, there are grey lines close to zero across the genome, suggesting that 
something else might fit.

![Figure 7: This sample does not match the expected parents, but it does look like it matches something else.](output/ibdpainting/2_H4_plot_ibd.png)

A quick look at the first rows of the score file for this sample suggests that the
combination of 175 and 6913 are likely candidates - the IBD score for this pair
is an order of magnitude lower than any other pair.
```
parent1 parent2 min_IBD
175     6913    0.013578783141498326
175     9743    0.12009568344735648
175     9852    0.12318103871182563
175     9581    0.1238256307101355
175     5210    0.1261425921515468
175     9972    0.12651983569807917
175     4840    0.12652044717407385
175     9821    0.12656418647709405
175     6924    0.12658849831752222
```
The other replicate of 9481x8307 shows the same pattern, and also matches 175 and 6913.
At the same time the two samples which should be derived from 175 and 6913 do not match those parents, but do match 9481 and 8307.

This looks like a simple case of 9481x8307 and 175x6913 being swapped.
Surprisingly, there is only one extra case where both parents have clearly been confused with something else (6012x997 F8 rep1 is probably 6021x9407), and one where a single parent has clearly been swapped for something else (9470x6107 F8 rep1
is probably 7002x6107).

#### Offspring is actually self-pollinated

Figure 8 shows what should be a cross between accessions 6108 and 6192.
However, it is in fact homozygous for the 6108 haplotype across the whole genome.
This is likely because the initial cross failed, and the line is in fact descended from a self-pollinated seed of 6108.
There is one additional case of this happening (6122x6252 replicate 2).

![Figure 8: Example of an F8 plant that is homozygous across the whole genome.](output/ibdpainting/1_B8_plot_ibd.png)

### Ambiguous cases

#### Missing haplotypes

Figure 9 shows a sample that matches an expected parent at most of the genome, but there are some regions of chromosomes 1, 2, 4 and 5 that don't match either parent, and don't really match anything else in the reference panel.

![Figure 9: This F8 plant matches the expected parent at most of the genome, but it is likely that the parts of the 6198 genome are still segregating.](output/ibdpainting/2_C11_plot_ibd.png)

The fact that *most* of the genome matches indicates that this is not a sample mix up.
There also don't appear to be any better candidates.
Here is the first ten lines of the `ibd_scores.csv` file:
```
parent1 parent2 min_IBD
6188    9381    0.07048910670002984
6197    9381    0.08793660411929506
6198    9381    0.08796521252504941
6188    9383    0.09783508908165152
6198    9383    0.1126233937330796
6197    9383    0.11262808158979185
6188    9382    0.13429055446514007
6191    9381    0.13632110527914162
6195    9381    0.14557521669338727
```
The expected parents (6198 and 9381) are there, and there are some other candidates with very similar names (6188 and 9383).
That might suggest that somebody had mixed up seed packets with very similar labels.
However, the plot shows that none of these are good matches.
Furthermore, the `min_IBD` scores are very similar for different candidate pairs; compare this table to the one [when labels were swapped](#labels-were-swapped), where the top candidate pairs had an order-of-magnitude lower IBD score than the rest.

This pattern is also unlikley to be residual heterozygosity left over from the cross; in the cases where there is heterozygosity (see above) there are usually hundreds of
heterozygous calls per megabase, but here there are an order of magnitude fewer.
You can also see this visually.
In genuinely heterozygous regions the lines for the expected parents are non-zero, but usually closer to zero than the point for the other accessions.
In this case the grey lines frequently *are* lower than those of the expected parents.

We think this pattern arises because of heterozygosity in the **parents** used to make the cross and reference panel, rather than arising from the cross itself.
The *A. thaliana* collection is a set of accessions sampled from natural populations, and although they are mostly self-pollinating, outcrossing does occur.
If an accession was sampled from the wild that was the ancestor of a cross some generations back then the heterozygous regions would still be segregating in the seed stock.
That would mean that the plant that was sampled for genotyping in the reference panel would not be the sample genotype as that used in the crosses, leading to partial mismatches.
This is even though both plants may have been homozygous themselves, and no human error has occured.
Consistent with this hypothesis, the replicate cross of the line shown in figure 9 shows mismatches in the same places (not plotted).

There are 45 lines that show pattern like this to a lesser or greater degree.
In these cases the material is correct, but the reference panel is not a good representation of the genotypes.
To get around this we will repeat the genotyping for these samples at high coverage to ensure SNPs are called accurately.

#### Parent appears incorrect

Sample 1_C4 matches expected parent 8283, but not 9471 (figure 10).
The stretches of mismatch are too long to be residual heterozygosity, and no other accessions in the reference panel match instead.
It is likely that some other accession was used in the cross, but that this accession was not in the reference panel.
Note that this is different from the case where [1137 and 1074 were missing](#expected-parents-are-not-in-the-reference-panel), because those crosses likely involve the correct parents.

![Figure 10: Sample that matches one parent, but the other parent appears to be missing.](output/ibdpainting/1_C4_plot_ibd.png)

There are 18 such cases where one parent appears to be missing, and one case where both parents seem to be missing.
Possible options here are to

1. Discard these samples
2. Treat them as valid and sequence at higher coverage

#### One parent *might* have been swapped

Sample 1_G10 matches expected parent 8256, but not 1367 (figure 11).
You can also see that there is at least one (multiple lines can hide each other) grey line close to zero in the regions where neither expected parent matches.
However, the grey lines are not really convincing matches.

![Figure 11: Sample that matches one but not both parents, and there is another candidate that partially matches the remaining haplotypes.](output/ibdpainting/1_G10_plot_ibd.png)

Similarly, a look at the IBD score file for this sample suggests that the missing parent might be 1061, but the evidence is not strong:
```
parent1 parent2 min_IBD
1061    8256    0.09085585092160826
1363    8256    0.13570305630994167
7522    8256    0.1444148987551863
1061    1363    0.1472248350711171
1066    8256    0.14796460020431212
6195    8256    0.1532090569124295
8256    9381    0.1634984095874054
5865    8256    0.16947001251860855
6092    8256    0.17021008221624542
```
Both the plot and the IBD score file for the other replicate of this cross show the same pattern.

So, there is evidence that the missing parent of this line is another accession in the reference panel, but the evidence is ambiguous.
An alternative explanation is that the true parent is really missing [as before](#parent-appears-incorrect), perhaps coupled with some residual or ancestral heterozygosity.
It is hard to know.
There are three additional crosses where both replicates look like one parent is missing and there is another candidate that *kinda* looks like it could be the parents, and two where only one replicate looks like there is an obvious other candidate (ten lines in total).

## Summary and action

Here is a summary of the outcomes of validation of 427 F8 lines:

1. 322 lines are good matches to expected parents. *No action required*.
2. 42 lines are likely matches to the expected parents, but the parent of the cross does not match the parent in the reference panel (because there is probably residual heterozygosity in the parental seed stock). *Sequence the F8s at higher coverage*.
3. 17 samples have no data. *Grow the plants again and repeat sequencing*.
4. 2 lines are self-pollinated. *Discard these lines*.
5. 8 lines have one or both parents that are likely correct, but the parents (1137 and 1074) are not in the reference panel. *Sequence 1137 and 1074 at high coverage*.
6. 6 lines do not match the expected parents, but are a close match to another accession in the reference panel. *Swap the labels*.
7. 11 lines do match one or both parents, and there is an ambiguous match to another accession in the reference panel. *Either resequence or discard*.
8. 22 lines do not match one or both parents. The true parent is likely not in the reference panel. *Sequence these at higher coverage*.