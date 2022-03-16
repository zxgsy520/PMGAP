#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import argparse
import logging

from pmgap.config import *
from thirdparty.dagflow import Task, DAG, do_dag
from pmgap.common import check_path, check_paths, cd, mkdir, get_version

LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "1.1.0"
__all__ = []

SOFTWARE_VERSION = {
    "abricate": {
        "GETVER": "%s/abricate --version 2>&1| grep 'abricate'" % ABRICATE_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.9.8"
    },
}


def create_abricate_task(gene, prefix, database, threads, job_type,
                         work_dir, out_dir):

    task = Task(
        id="%s_%s" % (database, prefix),
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={abricate}:$PATH
abricate {gene} --threads {threads} \\
  --mincov 50 --minid 75 \\
  --db {database}> {prefix}.{database}_old.tsv
{script}/abricate_summary.py {prefix}.{database}_old.tsv \\
  -o {prefix}.{database}_summary.tsv >{prefix}.{database}.tsv
cp {prefix}.{database}.tsv {prefix}.{database}_summary.tsv {out_dir}
""".format(id=id,
        abricate=ABRICATE_BIN,
        script=SCRIPTS,
        gene=gene,
        prefix=prefix,
        database=database,
        threads=threads,
        out_dir=out_dir)
    )

    return task, os.path.join(work_dir, '%s.%s_summary.tsv' % (prefix, database))


def create_abricate_tasks(gene, prefix, database, threads,
                          job_type, work_dir, out_dir):

    option = {}
    tasks =[]
    abricates = []
    for i in database.strip().split(';'):
        option[i] = {
            "version": get_version(SOFTWARE_VERSION["abricate"]),
            "option": "abricate  --mincov 50 --minid 75  --db %s" % i
        }
        temp_work = mkdir(os.path.join(work_dir, i))
        task, abricate = create_abricate_task(
            gene=gene,
            prefix=prefix,
            database=i,
            threads=threads,
            job_type=job_type,
            work_dir=temp_work,
            out_dir=out_dir)

        tasks.append(task)
        abricates.append(abricate)

    return tasks, abricates, option


def run_abricate(gene, prefix, threads, job_type, concurrent, refresh,
                work_dir="", out_dir="", database="card;vfdb"):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    gene = check_path(gene)
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    dag = DAG("abricate")
    tasks, abricates, option = create_abricate_tasks(
        gene=gene,
        prefix=prefix,
        database=database,
        threads=threads,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    options["software"].update(option)
    for i in tasks:
        dag.add_task(i)

    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    return options


def advanal_hlep_args(parser):

    parser.add_argument("gene", metavar='FILE', type=str, required=True,
        help="Input gene sequence file(fasta).")
    parser.add_argument("-p", "--prefix", metavar="STR", type=str, default="out",
        help="output file prefix.")
    parser.add_argument("-t", "--thread", metavar='INT', type=int, default=4,
        help="Set the running thread, default=4")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", default=".",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", default=".",
        help="Output directory (default: current directory)")

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
    run_abricate.py Annotation of drug-resistant genotoxicity factors.
attention:
    run_abricate.py gene.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = advanal_hlep_args(parser).parse_args()

    run_abricate(args.gene, args.prefix, args.threads, args.job_type, args.concurrent,
            args.refresh, args.work_dir, args.out_dir)


if __name__ == "__main__":

    main()
