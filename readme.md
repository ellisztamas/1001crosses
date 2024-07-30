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

Additional readme files in `01_data` give additional details as necessary.

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

#### Phenotyping experiment at 12°C

Tom Ellis set up a phenotyping experiment to collect data on rosette size (a 
proxy for growth), flowering time and seed size in the parents and F9s.

* 427 F9 families
* 219 parental lines
* Each line is replicated in three cohorts, with germination staggered by a few
weeks between cohorts.
* The experiment was meant to be at 10°C, but we realised the temperature sensor
in growth room 13 overreads, so we kept it at 12°C.
* We measured rosette size using the app easyLeafArea.
* Each pot label had a QR code encoding the position and ID of the plant. When
plants flowered I scanned them with an app that also recorded the date. Almudena
Morales also recorded some dates using a separate app.
* After flowering, plants were moved to greenhouse chamber 1 for seeds to mature.

Full details are documented on eLabJournal:
https://vbc.elabjournal.com/members/experiments/browser/#view=experiment&nodeID=253766&page=0&userID=0&status=0&column=created&order=DESC&search=

#### Flowering time (deprecated)

Joanna Gunis measured flowering time in the F8s.

We will probably not use these data because parents and F8s were not grown in the
same environment, there is on replication, and randomisation is not ideal.

#### External seed size data

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

#### Mamba environment

An environment file `environment.yml` is provided to handle (most) dependencies for running the scripts.
This does not include packages for R, as these do not seem to play well with conda/mamba.

I have used `micromamba` to manage dependencies, but conda would also work.
Assuming you have micromamba installed, install the environment with 
```
micromamba env create -f environment.yml
```

In the scripts in this repo this environment is loaded by `source setup.sh` 
which you will see at the top of most scripts.

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