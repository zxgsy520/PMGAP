#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def read_tsv(file, sep="\t"):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def stat_gff(gff, types):

    stat_dict = OrderedDict()

    for t in types:
        stat_dict[t] = [0, 0, '']

    for g in read_tsv(gff):
        if 'ncRNA_class=' in g[8]:
            _class = g[8].split('ncRNA_class=')[1].split(';')[0]
        else:
            _class = g[2]

        _type, _start, _end = g[2:5]

        if g[2] in stat_dict:
            stat_dict[_type][0] += 1
            stat_dict[_type][1] += abs(int(_start)-int(_end)) + 1
            stat_dict[_type][2] = _class

    print("#Type\tNumber\tLength (bp)")

    for k, v in stat_dict.items():
        print("{0}\t{1:,}\t{2:,}".format(v[2], v[0], v[1]))


def add_args(parser):

    parser.add_argument("gff", help="gff file")
    parser.add_argument("-t", "--types", nargs="+", metavar="TYPE", help="types to be stat")

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
    stat_gff(args.gff, args.types)


if __name__ == "__main__":
    main()

