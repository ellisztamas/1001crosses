# Linkage structure in *Arabidopsis thaliana*

A repo investigating the effect of one round of random mating on linkage
structure in *Arabidopsis thaliana*.

## Table of contents

1. [Experimental design](#experimental-design)
3. [Data](#data)
4. [Processing and analysis](#processing-and-analysis)
5. [Author information](#author-information)
6. [Acknowledgements](#acknowledgements)
7. [License](#license)

## Experimental design

## Data

### Design

Random pairs of 217 accessions from Sweden (that show strong population
structure) were crossed. Most of these crosses were reciprocal. The whole crossing
scheme was repeated to give two replicates of a bit more than 200 lines in 
duplicate. Each line was maintained by single-line descent for nine generations.

### Sequence data

#### Parental accessions

Sequence data for the parents are taken from a VCF file produced by Fernando
Rabanal from

> Brachi, Benjamin, et al. "Plant genetic effects on microbial hubs impact host fitness in repeated field trials."Proceedings of the National Academy of Sciences 119.30 (2022): e2201285119

This is a reanalysis of previously published data from the [1001 genomes project](https://1001genomes.org/)
updated with newer SNP calling methods.

#### F8 genotypes

429 lines were sequenced at the F8 generation by Pieter Clauw.
See NGS master list plates 2022-007 and 2022-008.

### Phenotypes

#### Flowering time

Joanna Gunis measured flowering time in the F8s.

Tom Ellis repeated this so that parents and offspring were grown in the same 
environment. Tom grew three replicates of all parents offspring (total = 637) in
Spring 2024 at 10°C.

#### Seed size

Tal Dahan sent seeds for the F8s and the parents (from the 2017 seed stock) to 
be measured by Boxeed in the Czech Republic. Since these cohorts are from
different maternal environments we also plan to send the F9s and parents from 
the flowering time experiment when these are mature.

## Processing and analysis

Scripts are separated into:
- `02_library`: reusable functions/applications.
- `03_processing`: data processing scripts to prepare raw data for further analysis
    and peform quality control.
- `04_output`: Deprecated.
- `05_results`: scripts to generate biological results.

All scripts are intended to be run from the project root folder.

Most of the scripts load a conda environment and set a **working directory** via
the script `setup.sh`.
This working directry is set by default to the path of the `scratch-cbe` drive
of the VBC CLIP cluster, which is optimised for large jobs with many temporary files.
This obviously won't work on a different machine, so change this path as
necessary. This can be done easily and only once by modifying `setup.sh`.

### Dependencies

#### Conda environment

A conda environment via the file `environment.yml` is provided with the necessary dependencies for running the scripts.
Install the environment with 
```
conda env create -f environment.yml
```

Near the start of bash scripts you will see a line that says something like
```
module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate 1001crosses
```
The first three lines of this load conda on the CLIP system.
Adjust as necessary for your machine.

#### R

This project uses R 4.2.1 and the tidyverse suite of packages.

## Author information

* Data collection:
    * Crosses by Fernando Rabanal and multipe assistants.
    * Maintenance of subsequent generations, phenotyping and tissue collection by Joanna Gunis.
    * Wet lab work by Viktoria Nyzhynska
* Analyses run by Tom Ellis based on previous work by Fernando Rabanal and Pieter Clauw.
* Principle investigator: Magnus Nordborg

## Acknoweldgements

We thank members of the Nordborg group over many years for constructive feedback about the project.