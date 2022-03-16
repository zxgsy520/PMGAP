#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

from Bio import SeqIO
from GffReader import open_gff


LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def _rev_comp(seq):

    r = ""

    transform_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    for i in seq:
        r += transform_dict[i]

    return r[::-1]


def get_gene(gff, genome, _type, gcode=11, translate=False):

    seqobj = SeqIO.index(genome, format="fasta")

    for g in open_gff(gff):
        seq = seqobj[g.seqid][g.start-1: g.end].seq
        
        if "N" in str(seq).upper():
            continue

        if g.strand == "-":
            seq = seq.reverse_complement()

        if g.type in _type:
            if g.type == "CDS" and translate:
                seq = seq.translate(table=gcode, cds=True)

            if 'ID' in g.attributes:
                _id = g.attributes["ID"]
            else:
                _id = g.attributes["Parent"]
            print(">%s [type=%s] [position=%s..%s,strand=%s]\n%s" % (_id, g.type, g.start, g.end, g.strand, seq))


def set_args():

    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("--gff", help="")
    args.add_argument("--genome", help="")
    args.add_argument("--type", nargs="+")
    args.add_argument("--gcode", type=int,)
    args.add_argument("--translate", action="store_true")

    return args.parse_args()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()
    get_gene(args.gff, args.genome, args.type, args.gcode, args.translate)


if __name__ == "__main__":
    main()

