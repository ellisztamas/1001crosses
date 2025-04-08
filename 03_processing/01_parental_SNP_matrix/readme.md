Scripts to pre-process SNPs for the parental lines.

SNP genotypes for 1163 inbred lines are in 
01_data/03_parental_genotypes/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf.gz
This includes all but two of the parental lines.
Accessions 1137 and 1074 are missing, and I sequenced them seprately.

These scripts calculate summary statistics on those SNPs, then filter for minor
allele count, missingness, quality, and depth.

This also creates a 'targets file' of variable SNP positions to call for cross
data later.