"""
Create a sample sheet for snpRarcher for the resequencing data.
The columns should be as follows (taken from the snpRarcher documentation):
    BioSample       The name of the sample. This will be the sample name in the final VCF
    LibraryName     LibraryID for sample, this can be the same or different than BioSample
    Run             The SRR for the sample, if applicable. If not, must be some unique value. It is often the lane number if samples are sequenced on multiple lanes.
    refGenome       Reference genome accession, if applicable.
    refPath         Path to local reference genome, if applicable.
    BioProject      If applicable. Otherwise any value is acceptable.
    fq1             Optional if no SRR value in Run. Path to read 1 for sample
    fq2             Optional if no SRR value in Run. Path to read 2 for sample
    SampleType      Optional. Triggers postproccesing module. Accepted values are ‘include’ or ‘exclude’

See the snpRarcher documentation for more information on these columns:
    https://snparcher.readthedocs.io/en/latest/setup.html

I'm omitting the long and lat columns as they are not relevant to this dataset.

Tom Ellis, 18th March 2024
"""
import pandas as pd
from glob import glob
import os
import re

import methlab as ml
print("Using methlab version " + ml.__version__)

# === Input ===

# Lists of paths to the fastq files separated by mate
path_to_fastq = "343561"
# path_to_fastq = "/scratch-cbe/users/thomas.ellis/crosses/04_resequencing/01_unzip_fastq/22YHLCLT3_4_R18523_20250303/demultiplexed/343561"
if os.path.exists(path_to_fastq):
    print(f"Using fastq files from {path_to_fastq}")
else:
    print(f"Error: could not find the path containing fastq files at {path_to_fastq}")
    print(f"Please check that you have unzipped the fastq files, and that you have the correct path.")
    print(f"The path is hard coded to a directory on the CLIP cluster, so you may need to change it.")
    print(f"Exiting...")
    exit()

mate1 = glob( path_to_fastq + "/*_R1_*.fastq.gz")
mate2 = glob( path_to_fastq + "/*_R2_*.fastq.gz")


# === Output ===

# CSV file matching the format required by snpRarcher
outfile="03_processing/04_resequencing/output/sample_sheet.csv"


# === Main ===

# Merge raw sample sheet  with paths to the fastq files.
sample_sheet = ml.align_fastq_with_plate_positions(mate1, mate2, raw_sample_sheet)

# Rename existing
sample_sheet = sample_sheet.rename(
    columns={
        'name': 'BioSample',
        'NGS sample ID': 'Run',
        'fq1' : 'fastq1',
        'fq2' : 'fastq2',
        })

# Add the other columns
sample_sheet['LibraryName'] = sample_sheet['BioSample']
sample_sheet['refGenome'] = 'TAIR10'
sample_sheet['refPath'] = '01_data/01_reference_genome/TAIR10_chr_all.fas'
sample_sheet['BioProject'] = '1001crosses'

# Reorder columns to what SNParcher expects
sample_sheet = sample_sheet[['BioSample', 'LibraryName', 'Run', 'refGenome', 'refPath', 'BioProject', 'fastq1', 'fastq2']]

# Write to file
sample_sheet.to_csv(outfile, index=False)
