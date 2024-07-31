"""
Create a list of new sample names to replace those in a VCF file.


Input:
    - CSV file giving plate number, well position and barcodes
    - A text file giving sample names in a VCF file to be replaced, with one sample
        per line (i.e. the output of `bcftools query -l vcf_file.vcf.gz`)
Output:
    A text file with the same number of rows as the input with new names for 
        each sample. New sample names are the plate number and well position of each
        sample.
"""

import pandas as pd
import os

# CSV file with sample IDs, plates, positions and barcodes.
plate_info = pd.read_csv(
    "01_data/02_F8_unaligned_bams/sequencing_plates_original.csv",
    dtype = str
    )
plate_info = plate_info.loc[~plate_info['name'].isna()]
# Single barcode entry. Notice the order is swapped.
plate_info['barcode'] = plate_info['barcode2'] + plate_info['barcode1'] 

# Text file with sample names to be replaced.
old_header = pd.read_csv(
    "/scratch-cbe/users/thomas.ellis/crosses/08_merge_VCF/header_names.txt",
    names = ['fullpath']
)
old_header['filename'] = [ os.path.basename(x) for x in old_header['fullpath'] ] # Remove paths
# Separate columns for NGS sample ID and barcode.
old_header[['NGS sample ID', 'barcode']] = old_header['filename'].str.split("_", expand = True)[[0,1]]

# Merge the two datasets
joined_df = pd.merge(
    old_header, plate_info , how = 'left', on = ['NGS sample ID', 'barcode']
)
# Column giving plate and position (e.g. 5_G6)
joined_df['new_name'] = joined_df['plate'] + "_" + joined_df['row'] + joined_df['col']

# Save to disk.
joined_df['new_name'].to_csv(
    "/scratch-cbe/users/thomas.ellis/crosses/08_merge_VCF/new_header_names.txt",
    index=False, header=False)