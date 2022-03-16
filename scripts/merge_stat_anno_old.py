#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        yield line.split('\t')


def read_connect(cstr):

    con_dict = {}
    
    for line in cstr.strip().split(';'):
        line = line.strip().split('=')
        con_dict[line[0]] = line[1]

    return con_dict


def read_gff(file):

    gff_dict = OrderedDict()

    for line in read_tsv(file):
        if line[2]!="gene":
            continue

        con_dict = read_connect(line[8])
        glen = int(line[4])-int(line[3])+1
        yield [con_dict["ID"], line[6], line[3], line[4], str(glen), line[0]]
        
        
def _combime(refseq, KEGG, COG, TIGRFAMs, Pfam, GO):

    stat_dict = {
        "Refseq":{},
        "Pfam":{},
        "TIGRFAMs":{},
        "COG":{},
        "KEGG":{},
        "GO":{},
        "name":{}
        }

    for record in sorted(read_tsv(refseq)):
        stat_dict["Refseq"][record[0]] = record[5].replace(",", "%2C").split(";")[0]

    for record in sorted(read_tsv(Pfam)):
        stat_dict["Pfam"][record[0]] = [record[2], record[5].replace(",", "%2C").split(";")[0]]

    for record in sorted(read_tsv(TIGRFAMs)):
        stat_dict["TIGRFAMs"][record[0]] = [record[2], record[5].replace(",", "%2C").split(";")[0]]
        if record[0] not in stat_dict["name"]:
            stat_dict["name"][record[0]] = record[4]
        else:
            if stat_dict["name"][record[0]] == "-":
                stat_dict["name"][record[0]] = record[4]
            else:
                stat_dict["name"][record[0]] += ",%s" % record[4]

    for record in sorted(read_tsv(COG)):
        stat_dict["COG"][record[0]] = [record[2], record[5].replace(",", "%2C").split(";")[0], record[3]]
        if record[0] not in stat_dict["name"]:
            stat_dict["name"][record[0]] = record[4]
        else:
            if stat_dict["name"][record[0]] == "-":
                stat_dict["name"][record[0]] = record[4]
            else:
                stat_dict["name"][record[0]] += ",%s" % record[4]

    for record in sorted(read_tsv(KEGG)):
        stat_dict["KEGG"][record[0]] = [record[2], record[5].replace(",", "%2C").split(";")[0], record[3]]
        if record[0] not in stat_dict["name"]:
            stat_dict["name"][record[0]] = record[4]
        else:
            if stat_dict["name"][record[0]] == "-":
                stat_dict["name"][record[0]] = record[4]
            else:
                stat_dict["name"][record[0]] += ",%s" % record[4]

    for record in sorted(read_tsv(GO)):

        goid = ""
        godes = ""

        for i in record[1:]:
            if i == "-":
                continue
            for j in i.split(";"):
                n, d = j.split(" - ")
                goid += "%s;" % n
                godes += "%s;" % d.replace(",", "%2C").replace("(", " - ").replace(")", "")
        
        stat_dict["GO"][record[0]] = [goid.strip(";"), godes.strip(";")]

    return stat_dict
    

def merge_anno(file, refseq, KEGG, COG, TIGRFAMs, Pfam, GO, out):
   
    output = open(out, "w")

    output.write("#Gene_Id\tStrand\tStart\tEnd\tGene_Length(bp)\tLocation\tGene_Name\tRefseq_Description\tPfam_Id\tPfam_Description\tTIGRFAMs_Id\tTIGRFAMs_Description\tCOG_Id\tCOG_Description\tCOG_Type\tKO_Id\tKO_Description\tPathway\tGO_Id\tGO_Description\n")

    stat_dict = _combime(refseq, KEGG, COG, TIGRFAMs, Pfam, GO)

    for line in read_gff(file):
        gstr = '\t'.join(line)

        if line[0] in stat_dict["name"]:
            gstr += "\t%s" % stat_dict["name"][line[0]]
        else:
            gstr += "\t-"

        if line[0] in stat_dict["Refseq"]:
            gstr += "\t%s" % stat_dict["Refseq"][line[0]]
        else:
            gstr += "\t-"
       
        if line[0] in stat_dict["Pfam"]:         
            gstr += "\t%s" % ('\t'.join(stat_dict["Pfam"][line[0]]))
        else:
            gstr += "\t%s" % ('\t'.join(["-", "-"]))

        if line[0] in stat_dict["TIGRFAMs"]:
            gstr += "\t%s" % ('\t'.join(stat_dict["TIGRFAMs"][line[0]]))
        else:
            gstr += "\t%s" % ('\t'.join(["-", "-"]))

        if line[0] in stat_dict["COG"]:
           gstr += "\t%s" % ('\t'.join(stat_dict["COG"][line[0]]))
        else:
           gstr += "\t%s" % ('\t'.join(["-", "-", "-"]))

        if line[0] in stat_dict["KEGG"]:
           gstr += "\t%s" % ('\t'.join(stat_dict["KEGG"][line[0]]))
        else:
           gstr += "\t%s" % ('\t'.join(["-", "-", "-"]))

        if line[0] in stat_dict["GO"]:
            gstr += "\t%s" % ('\t'.join(stat_dict["GO"][line[0]]))
        else:
            gstr += "\t%s" % ('\t'.join(["-", "-"]))

        output.write("%s\n" % gstr)
    output.close()


def add_help_args(parser):

    parser.add_argument("--gff", metavar='FILE', type=str, required=True,
        help="Gff file for input gene annotations.")
    parser.add_argument("--refseq", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the refseq database.")
    parser.add_argument("--KEGG", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the kegg database.")
    parser.add_argument("--COG", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the kegg database.")
    parser.add_argument("--TIGRFAMs", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the TIGRFAMs database.")
    parser.add_argument("--Pfam", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the Pfam database.")
    parser.add_argument("--GO", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the GO database.")
    parser.add_argument("-o", "--out", metavar='FLOAT', type=str, default="out.stat_anno.tsv",
        help="Output statistics table.")

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
    merge_stat_anno.py  Statistical summary of each database.

attention:
    merge_stat_anno.py -gff genomic.gff3 -rf refseq.tsv -kg KEGG.tsv -cog COG.tsv -o stat_anno.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    merge_anno(args.gff, args.refseq, args.KEGG, args.COG, args.TIGRFAMs, args.Pfam, args.GO, args.out)


if __name__ == "__main__":

    main()
