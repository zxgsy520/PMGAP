#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from collections import OrderedDict
from pmgap.config import *
from pmgap.parser import add_ngs_qc_args
from thirdparty.dagflow import DAG, Task, ParallelTask, do_dag
from pmgap.common import check_paths, mkdir, read_tsv, get_version


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "1131978210qq.com"
__version__ = "v1.2.0"


def merge_raw_data_task(name, r1, r2, job_type, work_dir, out_dir):

    suffix = ""
    tools = "cat"
    if r1[0].endswith(".gz") or r2[0].endswith(".gz"):
        suffix = ".gz"
        tools = "zcat"

    if len(r1)<=1:
        job_type = "local"
        run = """
ln -s {r1} {name}.raw.r1.fastq{suffix}
ln -s {r2} {name}.raw.r2.fastq{suffix}
""".format(r1=" ".join(r1), r2=" ".join(r2), name=name, suffix=suffix)
    else:
        suffix = ""
        run = """
{tools} {r1} >{name}.raw.r1.fastq
{tools} {r2} >{name}.raw.r2.fastq
""".format(tools=tools, r1=" ".join(r1), r2=" ".join(r2), name=name)

    task = Task(
        id="data_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{run}
""".format(run=run)
    )

    return task, os.path.join(work_dir, "%s.raw.r1.fastq%s" % (name, suffix)), os.path.join(work_dir, "%s.raw.r2.fastq%s" % (name, suffix))


def quality_control_task(r1, r2, name, trim, thread, job_type, work_dir, out_dir):

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
        id="ngs_qc_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={python}:{fastp}:{fastqc}:$PATH
fastp -i {r1} -I {r2} \
-o {name}.clean.r1.fastq -O {name}.clean.r2.fastq \
-w {thread} -n 0 -f {trim} -F {trim} -t {trim} -T {trim} --json {name}_fastp.json
fastqc {name}.clean.r1.fastq {name}.clean.r2.fastq -t {thread} --extract -o ./
python {scripts}/plot_fastqc.py -1 {name}.clean.r1_fastqc/fastqc_data.txt \
-2 {name}.clean.r2_fastqc/fastqc_data.txt --name {name}
python {scripts}/stat_fastp.py {name}_fastp.json --name {name} > {name}.qc.xls
cp {name}.qc.xls {name}.base_*.p* {out_dir}
""".format(scripts=SCRIPTS,
            python=PYTHON_BIN,
            fastp=FASTP_BIN,
            fastqc=FASTQC_BIN,
            thread=thread,
            name=name,
            r1=r1,
            r2=r2,
            trim=trim,
            out_dir=out_dir)
    )

    return task, os.path.join(work_dir, "%s.clean.r1.fastq" % name), os.path.join(work_dir, "%s.clean.r2.fastq" % name), option


def contamination_eval_task(name, r1, thread, job_type ,work_dir, out_dir):

    option = OrderedDict()
    option["blastn"] = {
        "version": get_version(SOFTWARE_VERSION["blastn"]),
        "option": "default"
    }

    task = Task(
        id="cont_eval_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={blast}:{python}:$PATH
python {scripts}/fq2fa.py {r1} -n 100000 >{name}.clean.r1.fa
blastn -query {name}.clean.r1.fa -db {dbase} \
-outfmt "6 std staxid sskingdom staxids" -max_target_seqs 5 -num_threads {thread} \
-out {name}.m6
python {scripts}/obtain_taxonomy.py -i {name}.m6 -t {taxonomy} -n {name}
python {scripts}/stat_taxonomy.py -i {name}.species_annotation.txt -rn 100000 -n {name}
python {scripts}/plot_stat_species.py -i {name}.stat_species.tsv -n {name}
cp {name}.species_classify.tsv {name}.stat_species.tsv {out_dir}
cp {name}.top10_species.tsv {name}.top10_species.png {name}.top10_species.pdf {out_dir}
cp {name}.species_annotation.txt {out_dir}
""".format(scripts=SCRIPTS,
            blast=BLAST_BIN,
            dbase=NT_DATABASE,
            python=PYTHON_BIN,
            taxonomy=TAXONOMY,
            name=name,
            r1=r1,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option


def ngs_qc_tasks(name, r1, r2, trim, thread, job_type, work_dir, out_dir):

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    merge_task, raw1, raw2 = merge_raw_data_task(
        name=name,
        r1=r1,
        r2=r2,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )

    qc_task, clean1, clean2, option = quality_control_task(
        r1=raw1,
        r2=raw2,
        name=name,
        trim=trim,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    qc_task.set_upstream(merge_task)
    options["software"].update(option)

    return merge_task, qc_task, options, clean1, clean2


def run_ngs_qc(r1, r2, name, trim, thread, job_type, concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    r1 = check_paths(r1)
    r2 = check_paths(r2)

    dag = DAG("ngs_qc")

    merge_task, qc_task, options, clean1, clean2 = ngs_qc_tasks(
        name=name,
        r1=r1,
        r2=r2,
        trim=trim,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )

    stat_qc = os.path.join(work_dir, "%s.qc.xls" % name)
    quality = os.path.join(work_dir, "%s.base_quality.png" % name)
    content = os.path.join(work_dir, "%s.base_content.png" % name)
    gc = os.path.join(work_dir, "%s.base_gc.png" % name)

    dag.add_task(merge_task)
    dag.add_task(qc_task)
    do_dag(dag, concurrent, refresh)

    with open(os.path.join(out_dir, "ngs_qc.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return clean1, clean2, stat_qc, quality, content, gc


def ngs_qc(args):

    run_ngs_qc(
        r1=args.read1,
        r2=args.read2,
        name=args.name,
        trim=args.trim,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )


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

    parser = add_ngs_qc_args(parser)
    args = parser.parse_args()
    ngs_qc(args)


if __name__ == "__main__":
    main()
