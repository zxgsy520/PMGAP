#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import collections

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("Reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    r = ""
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if (not line) or ("---" in line) or ("GENES" in line):
            continue
        if line.startswith(">"):
            if r:
                yield r.split("\n", 1)
            r = "%s\n" % line.strip(">")
            continue
        r += line.upper()

    if r:
        yield r.split("\n", 1)
    fp.close()


def tsv2gff(file, seqid, locus="NPGAP", gcode=11):

    n = 0
    r = {}
    for line in read_tsv(file):
        if len(line) <= 5 or line[0]=="No" or line[0]=="Stop":
            continue
        n += 1
        geneid = "%s_%05d" % (locus, n)
        print("{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={locus_tag}\n"
              "{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\tParent={locus_tag};transl_table={gcode}".format(
                seqid=seqid, source="vgas", type="CDS", start=line[1], end=line[2],
                score=line[5], strand=line[3], phase=".", locus_tag=geneid, gcode=gcode
                ))
        r[line[0]] = geneid

    return r


def seq2fasta(file, dictid, output):

    fo = open(output, "w")
    for seqid, seq in read_fasta(file):
        seqid = seqid.split()[2].split(":")[0]
        seqid = dictid[seqid]
        seqid = "{seqid} [gene={seqid}] [protein={seqid}]".format(seqid=seqid)
        fo.write(">%s\n%s\n" % (seqid, seq))

    fo.close()


def process_vgas(gff, gene, protein, genome, prefix="NG", locus="NG", gcode=11):

    if not locus:
        locus = prefix

    ids = []
    for seqid, seq in read_fasta(genome):
        ids.append(seqid.split()[0])
    if len(ids) >= 2:
        raise Exception("Too many genome sequences")

    r = tsv2gff(gff, ids[0], locus, gcode)
    seq2fasta(gene, r, "%s.gene.fasta" % prefix)
    seq2fasta(protein, r, "%s.protein.fasta" % prefix)

    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar='FILE', type=str,
        help="Input the gene location file annotated by vgas.")
    parser.add_argument("-g", "--gene", metavar="FILE", type=str, required=True,
        help="Input the gene sequence file annotated by vgas.")
    parser.add_argument("-p", "--protein", metavar="FILE", type=str, required=True,
        help="Input the protein sequence file annotated by vgas.")
    parser.add_argument("--genome", metavar="FILE", type=str, required=True,
        help="Input genome file.")
    parser.add_argument("--prefix", metavar="STR", type=str, default="NG",
        help="The prefix of the output file, default=NG.")
    parser.add_argument("-s", "--locus", metavar="STR", type=str, default="",
        help="Input the locus_tag of the gene, default=.")
    parser.add_argument("--gcode", metavar="INT", type=int, default=1,
        help="Gene codon number, default=1")

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
     process_vgas.py -- View multiple sequence alignment files

attetion:
    #muscle -align seq.fasta -output seq.aln.fasta
    add_look_alnfa_args.py seq.aln.fasta --ref EF641008.1 >seq.aln.xls

''')
    args = add_help_args(parser).parse_args()
    process_vgas(args.input, args.gene, args.protein, args.genome, args.prefix,
                 args.locus, args.gcode)


if __name__ == "__main__":
    main()
