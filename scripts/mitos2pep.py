#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("Reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def mitos2pep(file):

    seqids = []
    n = 0
    for line in read_tsv(file):
        if line[0].startswith(">"):
            gene = line[-1]
            #seqid = line[0].strip(";")
            if gene not in seqids:
                seqid = gene
                seqids.append(gene)
            else:
                n += 1
                seqid = "%s_%s" % (gene, n)
            print(">%s [gene=%s] [protein=%s]" % (seqid, gene, gene))
        else:
            #print(line[0].strip("*"))
            print(line[0])

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar='FILE', type=str,
        help="Input the protein sequence annotated by mitos2.")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    mitos2pep Format protein sequence

attention:
    mitos2pep protein.fasta >genomic.pep
    tbl2asn -i genomic.fasta -V b -s T

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    mitos2pep(args.input)


if __name__ == "__main__":

    main()
