#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from GffReader import open_gff


LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def read_gff(file):
    """

    [gene_gff, [CDS/RNA_gff]]

    :param file:
    :return:
    """
    pseqid = ""

    r = []

    for record in sorted(open_gff(file), key=lambda d: d.seqid):

        if record.seqid != pseqid:
            pseqid = record.seqid
            if r:
                yield r
                r = [record]
            print(">Feature %s" % pseqid)

        if record.type == "gene":
            if not r:
                r = [record]
                continue
            else:
                yield r
                r = [record]
        elif record.type in "CDS|rRNA|ncRNA|tRNA":
            r.append(record)
            continue
        else:
            yield [record]

    if r:
        yield r


def gff2tbl(file):

    for record in read_gff(file):

        if record[0].strand == "-":
            strand = "-"
        else:
            strand = "+"

        if len(record) == 1:

            if strand == "+":
                start, end = record[0].start, record[0].end
            else:
                start, end = record[0].end, record[0].start

            print("%s\t%s\t%s" % (start, end, record[0].type))

            for k in ["inference", "product"]:

                if k in record[0].attributes:
                    for i in record[0].attributes[k].split(","):
                        print("\t\t%s\t%s" % (k, i))
        else:
            if strand == "+":
                gstart, gend = record[0].start, record[0].end
                cstart, cend = record[0].start, record[0].end
            else:
                gstart, gend = record[0].end, record[0].start
                cstart, cend = record[0].end, record[0].start

            _attr = record[1].attributes
            locus_tag = _attr["Parent"]
            ctype = record[1].type

            if "gene" in _attr:
                gene = "\n\t\tgene\t%s" % _attr["gene"]
            else:
                gene = ""

            print("""\
{gstart}\t{gend}\tgene{gene}
\t\tlocus_tag\t{locus_tag}
{cstart}\t{cend}\t{ctype}{gene}
\t\tlocus_tag\t{locus_tag}""".format(**locals()))

            for k in ["EC_number", "inference", "note", "product"]:

                if k in record[1].attributes:
                    for i in record[1].attributes[k].split(","):
                        print("\t\t%s\t%s" % (k, i.replace("%2C", ",")))


def set_args():
    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("file", help="")

    return args.parse_args()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()

    gff2tbl(args.file)


if __name__ == "__main__":
    main()

