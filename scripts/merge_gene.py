#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn, invicoun@foxmail.com"
__all__ = []


def read_gff(file):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        r.append(line.split("\t"))

    return r


def merge_gene(files, locus="NPGAP", gcode=11):

    r = []
    for file in files:
        r += read_gff(file)

    n = 1
    is_n = 1

    print("##gff-version 3")
    for record in sorted(r, key=lambda d: (d[0], int(d[3]))):

        if record[1] == "PhiSpy" and record[2] != "prophage_region":
            continue
        elif record[2] == "CDS":
            if "partial=00" not in record[8]:
                continue

            record[1] = record[1].replace("_v", ":")
            print("{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={locus_tag}\n"
                  "{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\tParent={locus_tag};transl_table={gcode}".format(
                    seqid=record[0], source=record[1], type=record[2], start=record[3], end=record[4],
                    score=record[5], strand=record[6], phase=record[7], locus_tag="%s_%05d" % (locus, n), gcode=gcode
                    ))
        elif record[2] in "rRNA|tRNA":
            print("{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={locus_tag}\n"
                  "{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\tParent={locus_tag};{attr}".format(
                    seqid=record[0], source=record[1], type=record[2], start=record[3], end=record[4],
                    score=record[5], strand=record[6], phase=record[7], locus_tag="%s_%05d" % (locus, n), attr=record[8]
                    ))
        elif record[2] in "ncRNA|regulatory":
            continue
        elif record[2] in "genomic_island":
            print("{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\tID={id}".format(
                seqid=record[0], source=record[1], type=record[2], start=record[3], end=record[4],
                score=record[5], strand=record[6], phase=record[7], id="island_%s" % is_n
            ))

            is_n += 1
            continue
        else:
            print("\t".join(record))
            continue

        n += 1

    return 0


def set_args():
    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
name:
　　merge_gene　Combine gff files with genome structure annotations
attention:
　　merge_gene.py *.gff >genome.gff
version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("input", nargs="+", metavar="FILE", type=str,
        help="Input gff file for genome structure prediction.")
    args.add_argument("-s", "--locus", metavar="STR", type=str, default="NPGAP",
        help="Input the locus_tag of the gene, default=NPGAP.")
    args.add_argument("-g", "--gcode", metavar="INT", type=int, default=11,
        help="Gene codon number, default=11")

    return args.parse_args()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()

    merge_gene(args.input, args.locus, args.gcode)


if __name__ == "__main__":
    main()
