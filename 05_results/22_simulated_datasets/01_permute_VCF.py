"""
Permute genotypes at each marker in a VCF file

Randomise genotypes at each row in a VCF file to completely remove linkage
disequilibrium, even at short scales. 
This uses cyvcf2 (based on Cypthon) to parse the VCF rapidly.

Tom Ellis, 21st July 2025 using o3 High
"""

import argparse
import random
from cyvcf2 import VCF, Writer


def permute_vcf(in_vcf: str, out_vcf: str, seed: int | None = None) -> None:
    """Shuffle genotype columns for each variant independently."""
    if seed is not None:
        random.seed(seed)

    vcf_in = VCF(in_vcf)
    writer = Writer(out_vcf, vcf_in)  # header copied verbatim

    n_samples = len(vcf_in.samples)

    for var in vcf_in:                     # iterate over every locus
        gts = var.genotypes                # list[tuple]: (a1, a2, phased?)
        idx = list(range(n_samples))
        random.shuffle(idx)                # fresh permutation for THIS locus
        var.genotypes = [gts[i] for i in idx]
        writer.write_record(var)

    writer.close()
    vcf_in.close()



def main() -> None:
    ap = argparse.ArgumentParser(
        description="Permute genotype order at each locus in a VCF file."
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Input VCF (.vcf or .vcf.gz)")
    ap.add_argument("-o", "--output", required=True,
                    help="Output VCF (.vcf or .vcf.gz)")
    ap.add_argument("--seed", type=int,
                    help="Optional RNG seed for reproducibility")
    args = ap.parse_args()

    permute_vcf(args.input, args.output, args.seed)



if __name__ == "__main__":
    main()