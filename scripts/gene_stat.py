#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from collections import OrderedDict

from Bio import SeqIO


LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn, invicoun@foxmail.com"
__all__ = []


def gene_stat(gff, fasta):

    def read_tsv(file):

        for line in open(file):
            line = line.strip()

            if not line or line.startswith("#"):
                continue
            yield line.split("\t")

    gene_dict = OrderedDict([
        ("tRNA", []),
        ("rRNA", {}),
        ("CDS", []),
        ("CRISPR", []),
        ("genomic_island", []),
        ("prophage_region", []),
    ])

    for i in read_tsv(gff):
        _type, _start, _end = i[2], int(i[3]), int(i[4])
        _length = _end - _start + 1

        if _type == "rRNA":
            _product = i[-1].split("product=")[-1].split(";")[0].split()[0]

            if _product not in gene_dict["rRNA"]:
                gene_dict["rRNA"][_product] = []

            gene_dict["rRNA"][_product].append(_length)
        elif _type == "misc_feature":
            if i[-1].startswith("note=CRISPR array"):
                gene_dict["CRISPR"].append(_length)
        else:
            if _type in gene_dict:
                gene_dict[_type].append(_length)

    genome_length = sum([len(i) for i in SeqIO.parse(fasta, "fasta")])

    print("#Type\t\tNumber\tLength (bp)\t% genome")

    for k, v in gene_dict.items():
        if k == "rRNA":
            for i, j in sorted(v.items()):
                print("rRNA\t{0:}\t{1:,}\t{2:,}\t{3:.2f}".format(i, len(j), sum(j), sum(j) * 100.0 / genome_length))
        else:
            print("{0:}\t\t{1:,}\t{2:,}\t{3:.2f}".format(k, len(v), sum(v), sum(v)*100.0/genome_length))


def add_args(parser):

    parser.add_argument("--fasta", help="")
    parser.add_argument("--gff", help="")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""


version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()
    gene_stat(args.gff, args.fasta)


if __name__ == "__main__":
    main()

