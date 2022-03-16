#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from collections import OrderedDict


LOG = logging.getLogger(__name__)

__version__ = "1.1.1"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def read_tsv(file):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split("\t")


def _combime(r, refseq, SwissProt, KEGG, COG, TIGRFAMs, Pfam, GO):
    """

    :param refseq:
    :param KEGG:
    :param COG:
    :param TIGRFAMs:
    :param Pfam:
    :return: {id: {"inference": [], "product": [], "EC_number": [], "gene": [], "note":[],  }}
    """
    stat = OrderedDict([
        ("COG", set()),
        ("KEGG", set()),
        ("GO", set()),
        ("Refseq", set()),
        ("Pfam", set()),
        ("TIGRFAMs", set()),
        ("SwissProt", set())
    ])

    for record in sorted(read_tsv(GO)):
        stat["GO"].add(record[0])

        for k, v in zip(["GO_process", "GO_component", "GO_function"], record[1:]):

            if v != "-":
                for i in v.split(";"):
                    r[record[0]]["note"].append("%s: %s" % (i.replace(",", "%2C").replace("(", " - ").replace(")", "") , k))

    for record in sorted(read_tsv(TIGRFAMs)):
        stat["TIGRFAMs"].add(record[0])

        if record[4] != "-":
            r[record[0]]["gene"].append(record[4])

        r[record[0]]["inference"].append("protein motif:%s" % record[2])
        product = record[5].replace(",", "%2C").split(";")[0]
        r[record[0]]["note"].append("TIGRFAMs: %s" % product)
        r[record[0]]["product"].append(product)

        if record[10] != "-":
            r[record[0]]["EC_number"] += record[10].lstrip("EC:").split()

    for record in sorted(read_tsv(Pfam)):
        stat["Pfam"].add(record[0])
        r[record[0]]["inference"].append("protein motif:%s" % record[2])
        product = record[5].replace(",", "%2C").split(";")[0]
        r[record[0]]["note"].append("Pfam: %s" % product)
        r[record[0]]["product"].append(product)

    for record in sorted(read_tsv(COG)):
        stat["COG"].add(record[0])
        if record[4] != "-":
            r[record[0]]["gene"].append(record[4])
        cogid = record[2].replace(";"," ,")
        r[record[0]]["inference"].append("protein motif:%s" % cogid)
        product = record[5].replace(",", "%2C").split(";")[0]
        r[record[0]]["note"].append("COG: %s" % product)
        r[record[0]]["product"].append(product)

    for record in sorted(read_tsv(refseq)):
        stat["Refseq"].add(record[0])
        r[record[0]]["inference"].append("similar to AA sequence:%s" % record[2])
        product = record[5].replace(",", "%2C").split(";")[0]
        r[record[0]]["note"].append("Refseq: %s %s" % (record[1], product))
        r[record[0]]["product"].append(product)
        if record[4] != "-":
            r[record[0]]["gene"].append(record[4])

    for record in sorted(read_tsv(SwissProt)):
        stat["SwissProt"].add(record[0])
        #r[record[0]]["inference"].append("protein motif:%s" % record[2])
        product = record[5]
        #r[record[0]]["note"].append("SwissProt: %s" % product)
        r[record[0]]["product"].append(product)
        if record[4] != "-":
            r[record[0]]["gene"].append(record[4])
        
    for record in sorted(read_tsv(KEGG)):
        stat["KEGG"].add(record[0])
        if record[0] not in r:
            continue
        if record[4] != "-":
            r[record[0]]["gene"].append(record[4])
        r[record[0]]["inference"].append("similar to AA sequence:KEGG:%s" % (record[1]))
        product = record[5].replace(",", "%2C").split(";")[0]
        r[record[0]]["note"].append("KEGG: %s %s" % (record[1], product))
        r[record[0]]["product"].append(product)

        if record[10] != "-":
            r[record[0]]["EC_number"] += record[10].lstrip("EC:").split()

    return r, stat


def read_gff3(file):

    cds_dict = {}

    for record in read_tsv(file):
        _type = record[2]

        if _type == "CDS":
            _id = record[8].split(";")[0].split("=")[-1]
           
            cds_dict[_id] = {"inference": ["ab initio prediction:%s" % record[1]],
                             "product": [], "EC_number": [], "gene": [],
                             "note": []}

    return cds_dict


def confirm_product(annotate):

    product = "hypothetical protein"

    if not annotate["product"]:
        product = "hypothetical protein"
    else:
        for i in annotate["product"]:
            i = i.strip()
            if len(i) >= 1:
                 product = i
                 break

    return product
                     

def update_annotation(file, refseq, SwissProt, KEGG, COG, TIGRFAMs, Pfam, GO, prefix):

    cds_dict = read_gff3(file)

    cds_dict, stat_dict = _combime(cds_dict, refseq, SwissProt, KEGG, COG, TIGRFAMs, Pfam, GO)

    cds_num = 0
    print("##gff-version 3")
    for record in read_tsv(file):

        if record[2] == "CDS":

            cds_num += 1
            _id = record[8].split(";")[0].split("=")[-1]
            cds_dict[_id]["product"] = [confirm_product(cds_dict[_id])]

            for k in ["gene", "EC_number", "inference", "note", "product"]:

                if cds_dict[_id][k]:
                    if k == "gene":
                        record[-1] += ";%s=%s" % (k, cds_dict[_id][k][0])
                    elif k == "EC_number":
                        record[-1] += ";%s=%s" % (k, ",".join(set(cds_dict[_id][k])))
                    else:
                        record[-1] += ";%s=%s" % (k, ",".join(cds_dict[_id][k]))

        print("\t".join(record))

    fh = open("%s.function_summary.tsv" % prefix, "w")

    _all = set()
    _one = set()

    for k, v in stat_dict.items():
        fh.write("%s\t%s\t%.2f\n" % (k, len(v), len(v)*100.0/cds_num))

        if not _all:
            _all = v
        _all = _all & v

        _one = _one | v

    fh.write("""\
all databases\t%s\t%.2f
at least one databases\t%s\t%.2f
Overall\t%s\t100
""" % (len(_all), len(_all)*100.0/cds_num, len(_one), len(_one)*100.0/cds_num, cds_num))

    fh.close()


def add_help_args(parser):

    parser.add_argument("gff", metavar='FILE', type=str,
        help="Gff file for input gene annotations.")
    parser.add_argument("-sp", "--swissprot", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the SwissProt database.")
    parser.add_argument("-rf", "--refseq", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the refseq database.")
    parser.add_argument("-kg", "--KEGG", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the kegg database.")
    parser.add_argument("-cog", "--COG", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the kegg database.")
    parser.add_argument("-tf", "--TIGRFAMs", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the TIGRFAMs database.")
    parser.add_argument("-pf", "--Pfam", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the Pfam database.")
    parser.add_argument("-go", "--GO", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the GO database.")
    parser.add_argument("-p", "--prefix", metavar='FILE', type=str, default="out",
        help="output file prefix, default=out.")

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
    update_annotation.py  Statistical summary of each database.

attention:
    update_annotation.py genomic.gff3 -rf refseq.tsv -kg KEGG.tsv -cog COG.tsv >stat_anno.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    update_annotation(args.gff, args.refseq, args.swissprot, args.KEGG,
                      args.COG, args.TIGRFAMs, args.Pfam, args.GO, args.prefix
    )


if __name__ == "__main__":
    main()

