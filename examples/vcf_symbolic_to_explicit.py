#!/usr/bin/env python3

import argparse, pysam, sys
import logging
from enum import Enum


def parse_args():
    parser = argparse.ArgumentParser(description="VariantGrid API client")
    parser.add_argument('--fasta', required=True, help='Reference fasta file')
    parser.add_argument('--out', default='-', help='Output file (default: stdout)')
    parser.add_argument('--max-size', type=int, default=1000, help='Maximum variant size (default: 1000)')
    parser.add_argument('--genome-build', required=False, help='Genome build (valid for bioutils eg "GRCh37", "GRCh38"')
    parser.add_argument('input_vcf_filename')
    return parser.parse_args()


class FastaContigType(Enum):
    CHR_PREFIX = 1
    NO_CHR = 2
    ACCESSION = 3


class FastaChromWrapper:
    def __init__(self, fasta_filename, genome_build=None):
        """ Handles chrom conversion between VCF and fasta reference """
        self.fasta = pysam.FastaFile(fasta_filename)
        if "chr1" in self.fasta.references:
            self.fasta_contig_type = FastaContigType.CHR_PREFIX
        elif "1" in self.fasta.references:
            self.fasta_contig_type = FastaContigType.NO_CHR
        elif genome_build is None:
            raise ValueError(f"Fasta: '{fasta_filename}' doesn't have 'chr1' or '1' - you need to specify 'genome_build'")

        if genome_build:
            from bioutils.assemblies import make_name_ac_map
            self.name_ac_map = make_name_ac_map(genome_build)
            contig = self.name_ac_map["1"]
            if contig not in self.fasta.references:
                raise ValueError(f"Fasta: '{fasta_filename}' doesn't have '1' (mapped to '{contig}') from bioutils {genome_build=}")
            self.fasta_contig_type = FastaContigType.ACCESSION
        else:
            self.name_ac_map = None

    def _add_chr(self, chrom):
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        return chrom

    def _strip_chr(self, chrom):
        return chrom.replace("chr", "")

    def fetch(self, chrom, start, end):
        """ Handles all the chrom conversion """

        if self.fasta_contig_type == FastaContigType.CHR_PREFIX:
            fasta_chrom = self._strip_chr(chrom)
        elif self.fasta_contig_type == FastaContigType.NO_CHR:
            fasta_chrom = self._strip_chr(chrom)
        elif self.fasta_contig_type == FastaContigType.ACCESSION:
            lookup_chrom = self._strip_chr(chrom)
            fasta_chrom = self.name_ac_map[lookup_chrom]
        else:
            raise ValueError(f"Unknown fasta contig type: {self.fasta_contig_type}")
        return self.fasta.fetch(fasta_chrom, start, end)


if __name__ == "__main__":
    args = parse_args()
    fasta = FastaChromWrapper(args.fasta, args.genome_build)
    variant_file = pysam.VariantFile(args.input_vcf_filename)

    header = variant_file.header # .copy()
    if "LowUniqueAlignments" not in header.filters:  # hdr.filters is dict‑like
        header.add_line('##FILTER=<ID=LowUniqueAlignments,Description="Imported from TSO500">')

    if "GT" not in header.formats:
        # logging.error("No GT in header")
        header.add_meta('FORMAT',
                        items=[('ID', "GT"), ('Number', "."), ('Type', "String"), ('Description', "Genotype")])

    with pysam.VariantFile(args.out, mode="w", header=header) as vcf_out:
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
