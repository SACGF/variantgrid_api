#!/usr/bin/env python3

import argparse, pysam, sys
import logging


def parse_args():
    parser = argparse.ArgumentParser(description="VariantGrid API client")
    parser.add_argument('--fasta', required=True, help='Reference fasta file')
    parser.add_argument('--out', default='-', help='Output file (default: stdout)')
    parser.add_argument('--max-size', type=int, default=1000, help='Maximum variant size (default: 1000)')
    parser.add_argument('input_vcf_filename')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    fasta = pysam.FastaFile(args.fasta)
    variant_file = pysam.VariantFile(args.input_vcf_filename)

    header = variant_file.header # .copy()
    if "LowUniqueAlignments" not in header.filters:  # hdr.filters is dict‑like
        header.add_line('##FILTER=<ID=LowUniqueAlignments,Description="Imported from TSO500">')

    if "GT" not in header.formats:
        # logging.error("No GT in header")
        header.add_meta('FORMAT',
                        items=[('ID', "GT"), ('Number', "."), ('Type', "String"), ('Description', "Genotype")])

    with (pysam.VariantFile(args.out, mode="w", header=header) as vcf_out):
        # logging.error(f"HEADER.formats: {list(vcf_out.header.formats.items())}")

        for rec in variant_file:
            if rec.alts and rec.alts[0] in ("<DEL>", "<DUP>"):
                seq = fasta.fetch(rec.chrom, rec.pos - 1, rec.stop).upper()
                if rec.alts[0] == "<DEL>":
                    rec.ref = seq
                    rec.alts = (seq[0],)
                else:  # <DUP>
                    rec.ref = seq[0]
                    rec.alts = (seq,)

            if "GT" not in rec.format:
                # logging.error(f"No GT in format: {rec.format}")
                # logging.error(f"{dict(rec.format)}")
                for i, sample in enumerate(rec.samples):
                    # logging.error(f"{i}, {sample}")
                    fmt = rec.samples[i]
                    # logging.error(f"fmt: {type(fmt)}, {dir(fmt)}")
                    fmt["GT"] = tuple(None for _ in rec.alts)

            # logging.error(f"FORMAT: {dict(rec.format)}")

            try:
                vcf_out.write(rec)
            except OSError as e:
                logging.error("Bad record at %s:%s %s>%s END=%s",
                              rec.chrom, rec.pos, rec.ref, rec.alts, rec.info.get("END"))
                raise
