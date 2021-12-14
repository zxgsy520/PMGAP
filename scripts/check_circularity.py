#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
from Bio import SeqIO

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def do_blast(fasta, blastn):

    if not os.path.exists(blastn):
        raise Exception("blastn not found in %r" % blastn)

    cmd = "%s -task megablast -subject %s -query %s -evalue 1e-100 -gapopen 2 -gapextend 2 -outfmt 7" % (blastn, fasta, fasta)

    result = os.popen(cmd)

    for n, line in enumerate(result):
        line = line.strip()

        if line.startswith("#") or not line:
            continue

        lines = line.split("\t")

        if lines[0] == lines[1]:
            yield lines


def filter_blast_result(records, min_overlap_len, alt_overlap_ratio):

    r = {}

    for record in records:
        qseqid, sseqid, ident, aln_len, _, _, qstart, qend, sstart, send, _, _ = record
        ident, aln_len, qstart, qend, sstart, send = \
            float(ident), int(aln_len), int(qstart), int(qend), int(sstart), int(send)

        if qstart == sstart and qend == send:
            continue

        if aln_len < min_overlap_len:
            continue

        if ident < alt_overlap_ratio*100:
            continue

        if qseqid not in r:
            r[qseqid] = []
        r[qseqid].append([qseqid, ident, aln_len, qstart, qend, sstart, send])

    return r


def check_circularity(fasta, blastn, min_overlap_len, min_overlap_ratio, alt_overlap_ratio, max_overhang_len):

    blast_dict = filter_blast_result(
        do_blast(fasta, blastn),
        min_overlap_len=min_overlap_len,
        alt_overlap_ratio=alt_overlap_ratio
    )

    for record in SeqIO.parse(fasta, "fasta"):
        name = record.id

        if name not in blast_dict:
            record.description = "[topology=linear]"
            print(record.format("fasta"))
            continue

        circular = False

        for i in blast_dict[name]:
            qseqid, ident, aln_len, qstart, qend, sstart, send = i

            if ident < min_overlap_ratio*100:
                if not circular:
                    sys.stderr.write("%s Possibly circular. [%s-%s] matches [%s-%s] with %s%% identity.\n" % (
                        name, qstart, qend, sstart, send, ident))
            else:
                if qstart <= max_overhang_len and len(record) - send <= max_overhang_len:
                    sys.stderr.write("%s circular. [%s-%s] matches [%s-%s] with %s%% identity.\n" %
                                     (name, qstart, qend, sstart, send, ident))
                    circular = True
                    record.description = "[topology=circular] [completeness=complete]"
                    print(record[qstart-1:sstart-1].format("fasta"))
                    break

        if not circular:
            record.description = "[topology=linear]"
            print(record.format("fasta"))


def add_args(parser):

    parser.add_argument("fasta", help="")
    parser.add_argument("--min_overlap_len", type=int, default=1000, help="")
    parser.add_argument("--min_overlap_ratio", type=float, default=0.95, help="")
    parser.add_argument("--alt_overlap_ratio", type=float, default=0.95, help="")
    parser.add_argument("--max_overhang_len", type=float, default=100, help="")
    parser.add_argument("--blastn", default="/export/personal/software/software/blast/2.10.0/bin/blastn", help="")

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
    check_circularity(
        fasta=args.fasta,
        min_overlap_len=args.min_overlap_len,
        min_overlap_ratio=args.min_overlap_ratio,
        alt_overlap_ratio=args.alt_overlap_ratio,
        max_overhang_len=args.max_overhang_len,
        blastn=args.blastn
    )

if __name__ == "__main__":
    main()

