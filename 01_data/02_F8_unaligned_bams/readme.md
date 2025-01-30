Raw sequence data for the F8s.

See plates 2022-007 and 2022-008 on the [NGS master list](https://docs.google.com/spreadsheets/d/1XjO8zabj-1vlu-ex37MRnsnXB_c1U3_k-uXoeKaeGn0/edit#gid=1327165230).

Things to note:
- These are symlinks to the read-only files on /resources/ngs/nordborg/.
- The unzipped files will actually be fastq, not bam files.
- Plates 1-4 were run on cell H32TKDSX5 split across lanes 3 and 4 and need merging.
- Plate 5 was run on cell HMN2MDRX2.
- There are two sample sheets relating sample names to the barcodes used:
    - `sequencing_plates_original.csv` gives the names in the original sample
        sheet on the NGS master list. The exception is plate 5, which I took from
        Pieter's sample sheet because he realised that the rows of this plate
        are backwards.
    - `sequencing_plates_pieter.csv` gives Pieter's corrections after he ran SNPmatch on them.
- In plate 5, Pieter thinks G6 and H6 are the genotypes that were meant to be in
    B6 and A6 respectively. I added these manually to sequencing_plates_original.csv