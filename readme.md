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

### Sequence data

### Phenotypes

## Processing and analysis

Scripts are separated into:
- `02_library`: reusable functions/applications.
- `03_processing`: data processing scripts to prepare raw data for further analysis
    and peform quality control.
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

This Porject uses R 4.2.1.


tidyverse
googlesheets4

## Author information

* Data collection:
    * Crosses by Fernando Rabanal and multipe assistants.
    * Maintenance of subsequent generations, phenotyping and tissue collection by Joanna Gunis.
    * Wet lab work by Viktoria Nyzhynska
* Analyses run by Tom Ellis based on previous work by Fernando Rabanal and Pieter Clauw.
* Principle investigator: Magnus Nordborg

## Acknoweldgements

We thank members of the Nordborg group over many years for constructive feedback about the project.