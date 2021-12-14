#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from Bio import SeqIO
LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def circlator_proc(file, name):

    circular_seq = ""
    linear_seq = ""

    for record in SeqIO.parse(file, "fasta"):

        if "[topology=circular]" in record.description:
            circular_seq += record.format("fasta")
        else:
            linear_seq += record.format("fasta")

    if circular_seq:
        with open("%s.circular.fasta" % name, "w") as fh:
            fh.write(circular_seq)

    with open("%s.linear.fasta" % name, "w") as fh:
        fh.write(linear_seq)

    return 0


def add_args(parser):

    parser.add_argument("fasta", help="")
    parser.add_argument("name", help="")

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

    circlator_proc(args.fasta, args.name)


if __name__ == "__main__":
    main()

