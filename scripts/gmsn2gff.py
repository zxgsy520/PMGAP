#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_gff(file):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        line = line.split("\t")
        line[0] = line[0].split()[0]
        line[1] = line[1].split(".")[0]

        r.append(line)

    return r


def gmsn2gff(file, locus="NPGAP", gcode=1):

    r = read_gff(file)

    n = 1
    print("##gff-version 3")
    for record in sorted(r, key=lambda d: (d[0], int(d[3]))):
        if record[2] == "CDS":
            print("{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={locus_tag}\n"
                  "{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\tParent={locus_tag};transl_table={gcode}".format(
                    seqid=record[0], source=record[1], type=record[2], start=record[3], end=record[4],
                    score=record[5], strand=record[6], phase=record[7], locus_tag="%s_%05d" % (locus, n), gcode=gcode
                    ))
        else:
            print("\t".join(record))
            continue
        n += 1
    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input gff file for genome structure prediction.")
    parser.add_argument("-s", "--locus", metavar="STR", type=str, default="NPGAP",
        help="Input the locus_tag of the gene, default=NPGAP.")
    parser.add_argument("-g", "--gcode", metavar="INT", type=int, default=1,
        help="Gene codon number, default=1")

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
name:
　　gmsn2gff　Combine gff files with genome structure annotations
attention:
　　gmsn2gff.py *.gff >genome.gff
version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_help_args(parser)
    args = parser.parse_args()
    gmsn2gff(args.input, args.locus, args.gcode)


if __name__ == "__main__":
    main()
