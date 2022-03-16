#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn, 113178210@qq.com"
__all__ = []


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ""
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            if seq != "":
                yield seq.split('\n')
            seq = "%s\n" % line.strip(">")
        else:
            seq += line

    if seq != "":
        yield seq.split('\n')
    fp.close()


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    n = 0
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()
        if not line:
            continue
        n += 1
        if n ==1:
            seq.append(line.strip("@"))
            continue
        seq.append(line)
        if n ==4:
            n = 0
            seq = []
            yield seq
    fp.close()


def add_circular_edge(genome, minlen=10000):

    if genome.endswith(".fastq") or genome.endswith(".fq") or genome.endswith(".fastq.gz") or genome.endswith(".fq.gz"):
        fh = read_fastq(genome)
    elif genome.endswith(".fasta") or genome.endswith(".fa") or genome.endswith(".fasta.gz") or genome.endswith(".fa.gz"):
        fh = read_fasta(genome)
    else:
        raise Exception("%r file format error" % genome)

    r = ""
    for record in fh:
        seqid = record[0].split()[0]
        if "topology=circular" in record[0] or "circular=true" in record[0]:
            print(">%s [topology=circular] [completeness=complete]\n%s%s" % (seqid, record[1], record[1][:2500]))
        else:
            if len(record[1]) >= minlen:
                print(">%s [topology=linear]\n%s" % (seqid, record[1]))
            elif len(record[1]) >= 1000:
                r += ">%s\n%s\n" % (record[0], record[1])
            else:
                pass

    return r


def add_args(parser):

    parser.add_argument("fasta", help="")
    parser.add_argument("--minlen", metavar='INT', type=int, default=20000, help="")

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
    r = add_circular_edge(args.fasta, args.minlen)
    fo = open("unassembly.fasta", 'w')

    fo.write(r)
    fo.close()


if __name__ == "__main__":
    main()
