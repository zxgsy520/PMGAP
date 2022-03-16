#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []



def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fp.close()


def split_attr(attributes):

    r = {}
    contents = attributes.split(";")

    for content in contents:
        if not content:
            continue
        if "=" not in content:
            LOG.info("%r is not a good formated attribute: no tag!" % content)
            continue
        tag, value = content.split("=", 1)
        r[tag] = value

    return r


def to_string(attributes):

    attr = []

    for key, value in attributes.items():
        if key in "ID":
            attr.insert(0, '%s=%s' % (key, value))
        elif key in "Name":
            attr.insert(1, '%s=%s' % (key, value))
        elif key in "Parent":
            attr.insert(2, '%s=%s' % (key, value))
        else:
            attr.append('%s=%s' % (key, value))

    return ';'.join(attr)



def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">")
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line.split()[0])
            continue
        if len(seq) == 2:
            seq[1] += line
        else:
            seq.append(line)

    if len(seq) == 2:
        yield seq
    fp.close()



def read_geseq_gff(file):

    gene_dict = {}
    data = {}
    n = 0
    for line in read_tsv(file, "\t"):
        attrs = split_attr(line[-1])
        if "ID" not in attrs:
            continue
        geneid = "%s-%s" % (line[0], attrs["ID"])
        
        if "gene" == line[2]:
            if geneid in gene_dict:
                #raise Exception("Gene %s has a repeat" % geneid)
                continue
            gene_dict[geneid] = [int(line[3]), int(line[4]), line[6]]
            continue
        if line[2] not in ["CDS", "tRNA", "rRNA"]:
            continue

        if geneid not in gene_dict:
            LOG.info("%s" % attrs["ID"])
            continue
            gene_dict[geneid] = [int(line[3]), int(line[4]), line[6]]
        else:
            geneid = geneid.replace("%s-" %line[2].lower(), "gene-")
        if line[2] not in data:
            data[line[2]] = {}
        if geneid not in data[line[2]]:
            data[line[2]][geneid] = []

        data[line[2]][geneid].append([line[0], int(line[3]), int(line[4]), line[6]])

    return gene_dict, data


def complement(seq):

    cdict = {"A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    seq = list(seq.upper())
    nseq = ""
    for i in seq:
        nseq += cdict[i]

    return nseq


def reverse_complement(seq):

    seq = seq[::-1]

    return complement(seq)


def stat_geseq_gff(file, genome, prefix="OUT"):

    gene_dict, data = read_geseq_gff(file)
    genome_dict = {}
    genome_length = 0
    for seqid, seq in read_fasta(genome):
        genome_length += len(seq)
        genome_dict[seqid] = seq

    r = {}
    geneid = []
    for i in data:
        r[i] = []
        fo = open("%s.%s.fasta" % (prefix, i), "w")

        for j in data[i]:
            if i == "CDS":
                geneid.append(j)
            temp = []
            for seqid, start, end, position in data[i][j]:
                if start >= end:
                    start, end = end, start
                r[i].append(end-start+1)
                seq = genome_dict[seqid][start-1:end]
                if position == "-":
                    seq = reverse_complement(seq)
                temp.append(seq)
            if position == "-":
                temp = temp[::-1]
            fo.write(">%s\n%s\n" % (j, "".join(temp)))
        fo.close()

    print("#Type\tNumber\tTotal Length(bp)\tAverage Length(bp)\tProportion(%)")
    for i in r:
        print("{0}\t{1:,}\t{2:,}\t{3:,.2f}\t{4:.2f}".format(i, len(r[i]),
            sum(r[i]), sum(r[i])*1.0/len(r[i]), sum(r[i])*100.0/genome_length)
        )

    genes = []
    for i in geneid:
        start, end, position = gene_dict[i]
        if start >= end:
            start, end = end, start
        genes.append(end-start+1)

    print("{0}\t{1:,}\t{2:,}\t{3:,.2f}\t{4:.2f}".format("gene", len(genes),
        sum(genes), sum(genes)*1.0/len(genes), sum(genes)*100.0/genome_length)
    )


    return 0


def add_hlep_args(parser):

    parser.add_argument('gff', metavar='FILE', type=str,
        help='Input the gff file for geseq gene prediction.')
    parser.add_argument('-g', '--genome', metavar='FILE', type=str, required=True,
        help='Input genome file.')
    parser.add_argument('-p', '--prefix', metavar='STR', type=str, default="OUT",
        help='output file prefix.')

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
    stat_geseq_gff.py: Statistical CHLOE predicted gene structure
attention:
    stat_chloe_gff.py genomic.gff3 -g genome.fasta -p out >out.stat.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_geseq_gff(args.gff, args.genome, args.prefix)


if __name__ == "__main__":

    main()
