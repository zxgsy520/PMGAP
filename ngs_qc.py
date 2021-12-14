#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import random
import logging
import argparse

from collections import OrderedDict
from pmgap.config import *
from pmgap.parser import add_ngs_qc_args
from thirdparty.dagflow import DAG, Task, ParallelTask, do_dag
from pmgap.common import check_paths, mkdir, read_tsv, get_version


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.2.1"


def merge_data_task(prefix, reads1, reads2, job_type, work_dir, out_dir):

    suffix = ""
    tools = "cat"
    if reads1[0].endswith(".gz") or reads2[0].endswith(".gz"):
        suffix = ".gz"
        tools = "zcat"

    if len(reads1)<=1:
        job_type = "local"
        run = """
ln -s {read1} {prefix}.raw.r1.fastq{suffix}
ln -s {read2} {prefix}.raw.r2.fastq{suffix}
""".format(read1=" ".join(reads1), read2=" ".join(reads2),
           prefix=prefix, suffix=suffix)
    else:
        suffix = ""
        run = """
{tools} {read1} >{prefix}.raw.r1.fastq
{tools} {read2} >{prefix}.raw.r2.fastq
""".format(tools=tools, read1=" ".join(reads1),
           read2=" ".join(reads2), prefix=prefix)

    task = Task(
        id="merge_data",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
{run}
""".format(run=run)
    )

    raw_r1 = os.path.join(work_dir, "%s.raw.r1.fastq%s" % (prefix, suffix))
    raw_r2 = os.path.join(work_dir, "%s.raw.r2.fastq%s" % (prefix, suffix))

    return task, raw_r1, raw_r2


def quality_control_task(read1, read2, prefix, trim, thread,
                         job_type, work_dir, out_dir):

    option = OrderedDict()
    option["fastp"] = {
        "version": get_version(SOFTWARE_VERSION["fastp"]),
        "option": "-n 0 -f {trim} -F {trim} -t {trim} -T {trim}".format(trim=trim)
    }
    option["fastqc"] = {
        "version": get_version(SOFTWARE_VERSION["fastqc"]),
        "option": "--extract"
    }

    task = Task(
        id="ngs_qc",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={python}:{fastp}:{fastqc}:$PATH
fastp -i {read1} -I {read2} \
-o {prefix}.clean.r1.fastq -O {prefix}.clean.r2.fastq \
-w {thread} -n 0 -f {trim} -F {trim} -t {trim} -T {trim} --json {prefix}_fastp.json
fastqc {prefix}.clean.r1.fastq {prefix}.clean.r2.fastq -t {thread} --extract -o ./
python {scripts}/plot_fastqc.py -1 {prefix}.clean.r1_fastqc/fastqc_data.txt \
-2 {prefix}.clean.r2_fastqc/fastqc_data.txt --name {prefix}
python {scripts}/stat_fastp.py {prefix}_fastp.json --name {prefix} > {prefix}.qc.xls
cp {prefix}.qc.xls {prefix}.base_*.p* {out_dir}
""".format(scripts=SCRIPTS,
            python=PYTHON_BIN,
            fastp=FASTP_BIN,
            fastqc=FASTQC_BIN,
            thread=thread,
            prefix=prefix,
            read1=read1,
            read2=read2,
            trim=trim,
            out_dir=out_dir)
    )
    clean_r1 = os.path.join(work_dir, "%s.clean.r1.fastq" % prefix)
    clean_r2 = os.path.join(work_dir, "%s.clean.r2.fastq" % prefix)

    return task, clean_r1, clean_r2, option


def blast_read_task(prefix, read1, thread, job_type, work_dir, out_dir):

    option = OrderedDict()
    option["blastn"] = {
        "version": get_version(SOFTWARE_VERSION["blastn"]),
        "option": "default"
    }

    task = Task(
        id="blast",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={blast}:{python}:$PATH
python {scripts}/fq2fa.py {read1} -n 100000 >{prefix}.clean.r1.fa
blastn -query {prefix}.clean.r1.fa -db {dbase} \\
  -outfmt "6 std staxid sskingdom staxids" -max_target_seqs 5 -num_threads {thread} \\
  -out {prefix}.m6
python {scripts}/obtain_taxonomy.py -i {prefix}.m6 -t {taxonomy} -n {prefix}
python {scripts}/stat_taxonomy.py -i {prefix}.species_annotation.txt -rn 100000 -n {prefix}
python {scripts}/plot_stat_species.py {prefix}.stat_species.tsv -p {prefix} > {prefix}.top10_species.tsv
cp {prefix}.species_classify.tsv {prefix}.stat_species.tsv {out_dir}
cp {prefix}.top10_species.png {prefix}.top10_species.pdf {out_dir}
cp {prefix}.species_annotation.txt {prefix}.top10_species.tsv {out_dir}
""".format(scripts=SCRIPTS,
            blast=BLAST_BIN,
            dbase=NT_DATABASE,
            python=PYTHON_BIN,
            taxonomy=TAXONOMY,
            prefix=prefix,
            read1=read1,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option


def run_ngs_qc(reads1, reads2, prefix, trim, thread, job_type,
               concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    reads1 = check_paths(reads1)
    reads2 = check_paths(reads2)

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    dag = DAG("ngs_qc")

    merge_task, raw1, raw2 = merge_data_task(
        prefix=prefix,
        reads1=reads1,
        reads2=reads2,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(merge_task)

    qc_task, clean1, clean2, option = quality_control_task(
        read1=raw1,
        read2=raw2,
        prefix=prefix,
        trim=trim,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(qc_task)
    qc_task.set_upstream(merge_task)
    options["software"].update(option)

    x = random.randint(0, len(reads1)-1)
    blast_task, option = blast_read_task(
        prefix=prefix,
        read1=reads1[x],
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(blast_task)
    options["software"].update(option)
    do_dag(dag, concurrent, refresh)

    stat_qc = os.path.join(work_dir, "%s.qc.xls" % prefix)
    quality = os.path.join(work_dir, "%s.base_quality.png" % prefix)
    content = os.path.join(work_dir, "%s.base_content.png" % prefix)
    gc = os.path.join(work_dir, "%s.base_gc.png" % prefix)

    return clean1, clean2, stat_qc, quality, content, gc, options


def ngs_qc(args):

    clean1, clean2, stat_qc, quality, content, gc, options = run_ngs_qc(
        reads1=args.reads1,
        reads2=args.reads2,
        prefix=args.prefix,
        trim=args.trim,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )
    with open(os.path.join(args.out_dir, "ngs_qc.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return 0


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
describe:
    Second-generation data quality control process

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_ngs_qc_args(parser)
    args = parser.parse_args()
    ngs_qc(args)


if __name__ == "__main__":
    main()
