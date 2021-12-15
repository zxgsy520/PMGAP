#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    seq = ''
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            if seq:
                seq = seq.split('\n')
                yield seq[0], seq[1]
            seq = "%s\n" % line.strip(">")
        else:
            seq += line

    if seq:
        seq = seq.split('\n')
        yield seq[0], seq[1]

    fp.close()


def n50(lengths):

    sum_length = sum(lengths)
    accu_length = 0
    nl = 0

    for i in sorted(lengths, reverse=True):
        accu_length += i
        if accu_length >= sum_length*0.5:
            nl = i
            break

    return nl


def stat_genome(genome):

    lengths = []
    cirs = 0

    for seqid, seq in read_fasta(genome):
        lengths.append(len(seq))
        if "circular" in seqid:
            cirs += 1

    return sum(lengths), len(lengths), n50(lengths), max(lengths), min(lengths), cirs


def get_sample(file):

    sample = file.split("/")[-1]
    sample = sample.split(".")[0]

    return sample


def stat_genomes(genomes):

    print("#Sample\tTotal bases\tContig number\tContig N50\tLongest Contig\tShortest Contig\tCircular Contig")

    for genome in genomes:
        sample = get_sample(genome)
        tb, cn, cn50, lc, sc, cirs = stat_genome(genome)
        print("{0}\t{1:,}\t{2:,}\t{3:,}\t{4:,}\t{5:,}\t{6:,}".format(sample,
              tb, cn, cn50, lc, sc, cirs))

    return 0


def add_help_args(parser):

    parser.add_argument('genomes', nargs='+', metavar='FILE', type=str,
        help='Input the genome file (fasta).')

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
     stat_genomes -- Statistical genome length and base content

attention:
    stat_genomes *.genome.fasta >stat_genomes.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))
    args = add_help_args(parser).parse_args()

    stat_genomes(args.genomes)


if __name__ == "__main__":
    main()
