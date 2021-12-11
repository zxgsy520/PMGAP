#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import random
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "2.0.0"
__author__ = ("invicoun@foxmail.com",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fastq(file):
    '''Read fastq file'''
    if file.endswith(".fq.gz") or file.endswith(".fastq.gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        seq.append(line)
        if len(seq)>=4:
            seq[0] = seq[0].strip("@").split()[0]
            yield seq[0], seq[1]
            seq = []

    fp.close()


def cut_seq(seq, minlen=250, maxlen=500):

    seqlen = len(seq)

    if seqlen > maxlen:
        lenseq = random.randint(minlen, maxlen)
        end = random.randint(lenseq, seqlen)
        start = end-lenseq
        seq = seq[start:end]

    return seq


def fq2fa(file, number, cut):

    n = 0

    if number!='all':
        number=int(number)

    for seq_id,seq in read_fastq(file):
        if n==number:
             break
        if cut:
            seq = cut_seq(seq)

        print('>%s\n%s' % (seq_id,seq))
        n +=1


def add_args(parser):

    parser.add_argument('fastq', help='Input fastq file.')
    parser.add_argument('-n', '--number', default='all',
        help='Output reads number, default=all')
    parser.add_argument("--cut", action="store_true",
        help="Cutting sequence.")
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
    fq2fa.py Extract sequence

attention:
    fq2fa.py txt.fastq >txt.fasta
    fq2fa.py txt.fastq -n 10000 >txt.fasta
    fq2fa.py txt.fastq -n 10000 --cut >txt.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))
    args = add_args(parser).parse_args()

    fq2fa(args.fastq, args.number, args.cut)


if __name__ == "__main__":
    main()
