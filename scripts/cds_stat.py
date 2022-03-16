#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from collections import OrderedDict
from Bio import SeqIO, Seq


LOG = logging.getLogger(__name__)

__version__ = "0.1.1"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn, invicoun@foxmail.com"
__all__ = []


def codon_usage(file, name, gcode=11, min_length=300):

    r = []

    codon = OrderedDict([
        ("GGG", [Seq.Seq("GGG").translate(table=gcode), 0]),
        ("GGA", [Seq.Seq("GGA").translate(table=gcode), 0]),
        ("GGT", [Seq.Seq("GGT").translate(table=gcode), 0]),
        ("GGC", [Seq.Seq("GGC").translate(table=gcode), 0]),
        ("GAG", [Seq.Seq("GAG").translate(table=gcode), 0]),
        ("GAA", [Seq.Seq("GAA").translate(table=gcode), 0]),
        ("GAT", [Seq.Seq("GAT").translate(table=gcode), 0]),
        ("GAC", [Seq.Seq("GAC").translate(table=gcode), 0]),
        ("GTG", [Seq.Seq("GTG").translate(table=gcode), 0]),
        ("GTA", [Seq.Seq("GTA").translate(table=gcode), 0]),
        ("GTT", [Seq.Seq("GTT").translate(table=gcode), 0]),
        ("GTC", [Seq.Seq("GTC").translate(table=gcode), 0]),
        ("GCG", [Seq.Seq("GCG").translate(table=gcode), 0]),
        ("GCA", [Seq.Seq("GCA").translate(table=gcode), 0]),
        ("GCT", [Seq.Seq("GCT").translate(table=gcode), 0]),
        ("GCC", [Seq.Seq("GCC").translate(table=gcode), 0]),
        ("AGG", [Seq.Seq("AGG").translate(table=gcode), 0]),
        ("AGA", [Seq.Seq("AGA").translate(table=gcode), 0]),
        ("AGT", [Seq.Seq("AGT").translate(table=gcode), 0]),
        ("AGC", [Seq.Seq("AGC").translate(table=gcode), 0]),
        ("AAG", [Seq.Seq("AAG").translate(table=gcode), 0]),
        ("AAA", [Seq.Seq("AAA").translate(table=gcode), 0]),
        ("AAT", [Seq.Seq("AAT").translate(table=gcode), 0]),
        ("AAC", [Seq.Seq("AAC").translate(table=gcode), 0]),
        ("ATG", [Seq.Seq("ATG").translate(table=gcode), 0]),
        ("ATA", [Seq.Seq("ATA").translate(table=gcode), 0]),
        ("ATT", [Seq.Seq("ATT").translate(table=gcode), 0]),
        ("ATC", [Seq.Seq("ATC").translate(table=gcode), 0]),
        ("ACG", [Seq.Seq("ACG").translate(table=gcode), 0]),
        ("ACA", [Seq.Seq("ACA").translate(table=gcode), 0]),
        ("ACT", [Seq.Seq("ACT").translate(table=gcode), 0]),
        ("ACC", [Seq.Seq("ACC").translate(table=gcode), 0]),
        ("TGG", [Seq.Seq("TGG").translate(table=gcode), 0]),
        ("TGA", [Seq.Seq("TGA").translate(table=gcode), 0]),
        ("TGT", [Seq.Seq("TGT").translate(table=gcode), 0]),
        ("TGC", [Seq.Seq("TGC").translate(table=gcode), 0]),
        ("TAG", [Seq.Seq("TAG").translate(table=gcode), 0]),
        ("TAA", [Seq.Seq("TAA").translate(table=gcode), 0]),
        ("TAT", [Seq.Seq("TAT").translate(table=gcode), 0]),
        ("TAC", [Seq.Seq("TAC").translate(table=gcode), 0]),
        ("TTG", [Seq.Seq("TTG").translate(table=gcode), 0]),
        ("TTA", [Seq.Seq("TTA").translate(table=gcode), 0]),
        ("TTT", [Seq.Seq("TTT").translate(table=gcode), 0]),
        ("TTC", [Seq.Seq("TTC").translate(table=gcode), 0]),
        ("TCG", [Seq.Seq("TCG").translate(table=gcode), 0]),
        ("TCA", [Seq.Seq("TCA").translate(table=gcode), 0]),
        ("TCT", [Seq.Seq("TCT").translate(table=gcode), 0]),
        ("TCC", [Seq.Seq("TCC").translate(table=gcode), 0]),
        ("CGG", [Seq.Seq("CGG").translate(table=gcode), 0]),
        ("CGA", [Seq.Seq("CGA").translate(table=gcode), 0]),
        ("CGT", [Seq.Seq("CGT").translate(table=gcode), 0]),
        ("CGC", [Seq.Seq("CGC").translate(table=gcode), 0]),
        ("CAG", [Seq.Seq("CAG").translate(table=gcode), 0]),
        ("CAA", [Seq.Seq("CAA").translate(table=gcode), 0]),
        ("CAT", [Seq.Seq("CAT").translate(table=gcode), 0]),
        ("CAC", [Seq.Seq("CAC").translate(table=gcode), 0]),
        ("CTG", [Seq.Seq("CTG").translate(table=gcode), 0]),
        ("CTA", [Seq.Seq("CTA").translate(table=gcode), 0]),
        ("CTT", [Seq.Seq("CTT").translate(table=gcode), 0]),
        ("CTC", [Seq.Seq("CTC").translate(table=gcode), 0]),
        ("CCG", [Seq.Seq("CCG").translate(table=gcode), 0]),
        ("CCA", [Seq.Seq("CCA").translate(table=gcode), 0]),
        ("CCT", [Seq.Seq("CCT").translate(table=gcode), 0]),
        ("CCC", [Seq.Seq("CCC").translate(table=gcode), 0]),
        ]
    )

    aa_dict = {
        "A": "Ala",
        "R": "Arg",
        "N": "Asn",
        "D": "Asp",
        "C": "Cys",
        "Q": "Gln",
        "E": "Glu",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "L": "Leu",
        "K": "Lys",
        "M": "Met",
        "F": "Phe",
        "P": "Pro",
        "S": "Ser",
        "T": "Thr",
        "W": "Trp",
        "Y": "Tyr",
        "V": "Val",
        "*": "*"
    }

    for record in SeqIO.parse(file, "fasta"):

        if len(record) % 3 != 0:
            continue

        if "=CDS" in record.description:
            r.append(len(record)/3)

            if len(record) < min_length:
                continue

            for i in [record.seq[i:i+3] for i in range(0, len(record), 3)]:
                codon[i][1] += 1

    codon_sum = sum([i[1] for i in codon.values()])

    fh = open("%s.codon.tsv" % name, "w")
    fh.write("#AA\tCodon\tNumber\t/1000\n")

    for c in codon:
        aa, num = codon[c]
        fh.write("%s\t%s\t%s\t%.2f\n" % (aa_dict[aa], c, format(num, ","), num*1000.0/codon_sum))

    fh.close()

    return r


def plot_pep(lengths, window, name):

    x = [(i+0.5)*window for i in range(int(max(lengths)/window)+1)]
    y = [0 for _ in x]

    for i in lengths:
        y[int(i/window)] += 1

    import matplotlib
    matplotlib.use("Agg")    
    from matplotlib import pyplot as plt
    #plt.style.use('seaborn')
    fig, ax = plt.subplots(figsize=(6, 4.5), )
    plt.subplots_adjust(top=0.95, left=0.10, right=0.95, bottom=0.10)
    ax.bar(x, y, width=window, linewidth=0.5, edgecolor="white", )
    maxy = window*-0.1
    if maxy <= 0:
        maxy = 0
    plt.ylim([maxy, plt.ylim()[1]])
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.xlabel("Protein length (aa)", fontsize=10, weight="bold")
    plt.ylabel("Number", weight="bold", fontsize=10,)
    plt.savefig("%s.protein_length.pdf" % name)
    plt.savefig("%s.protein_length.png" % name, dpi=900)


def set_args():

    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("fasta", help="")
    args.add_argument("--name", type=str, required=True, help="")
    args.add_argument("--min_length", type=int, default=0, help="")
    args.add_argument("--gcode", type=int, default=11, help="")

    return args.parse_args()


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()
    lengths = codon_usage(args.fasta, args.name, args.gcode, args.min_length)
    plot_pep(lengths, 100, args.name)


if __name__ == "__main__":
    main()

