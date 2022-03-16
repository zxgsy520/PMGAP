#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from pmgap.config import *
from thirdparty.dagflow import DAG, Task, do_dag
from pmgap.common import check_path, mkdir, get_version
from pmgap.parser import add_get_mito_args

LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.2.0"


def get_ngs_reads_task(prefix, database, platform, read1, read2, job_type,
                       work_dir, out_dir, thread=4, gc=20, base="all", score=0):

    option = {}
    option["recycle"] = {
        "version": "v1.0.0",
        "option": "choose_ml_map --gc %s --base %s --score %s" % (gc, base, score)
    }

    task=Task(
        id='get_ngs',
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={minimap2}:$PATH
minimap2 -t {thread} -x sr {database} {read1} {read2} >{prefix}.ngs.paf
{bin}/choose_ml_map --input {prefix}.ngs.paf --platform {platform} --gc {gc} --base {base} --score {score} \\
  --read1 {read1} --read2 {read2} --prefix {prefix}
 #cp {prefix}.choose_r1.fq {prefix}.choose_r2.fq
""".format(minimap2=MINIMAP2_BIN,
            bin=BIN,
            database=database,
            prefix=prefix,
            platform=platform,
            read1=read1,
            read2=read2,
            gc=gc,
            base=base,
            score=score,
            thread=thread
        )
    )

    mito_r1 =  os.path.join(work_dir, '%s.choose_r1.fq' % prefix)
    mito_r2 =  os.path.join(work_dir, '%s.choose_r2.fq' % prefix)

    return task, mito_r1, mito_r2, option


def get_tgs_reads_task(prefix, database, reads, job_type, work_dir, out_dir,
                       thread=4, gc=20, base="all", score=5, sequencer=""):

    task=Task(
        id='get_tgs',
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={minimap2}:$PATH
minimap2 -t {thread} {x} {database} {reads} >{prefix}.tgs.paf
{bin}/choose_ml_map --input {prefix}.tgs.paf --platform illumina --gc {gc} --base {base} --score {score} \\
  --read1 {reads} --prefix {prefix}
mv {prefix}.choose_r1.fq {prefix}.choose.fq
#cp {prefix}.choose.fq {out_dir}
""".format(minimap2=MINIMAP2_BIN,
            bin=BIN,
            x=SEQUENCER[sequencer]["minimap2"].replace("-ax", "-x"),
            database=database,
            prefix=prefix,
            reads=reads,
            gc=gc,
            base=base,
            score=score,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, '%s.choose.fq' % prefix)


def run_get_mito(prefix, database, platform, read1, read2, sequencer, reads,
                 job_type, work_dir, out_dir, concurrent, refresh,
                 thread=4, gc=20, base="all", score=0):

    read1 = check_path(read1)
    read2 = check_path(read2)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    database = check_path(database)

    work_dict = {
        "ngs": "00_get_ngs",
        "tgs": "01_get_tgs"
    }

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    dag = DAG("mit_assembly")
    ngs_task,  mito_r1, mito_r2, option = get_ngs_reads_task(
        prefix=prefix,
        database=database,
        platform=platform,
        read1=read1,
        read2=read2,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["ngs"]),
        out_dir=out_dir,
        thread=thread,
        gc=gc,
        base=base,
        score=score
    )
    dag.add_task(ngs_task)
    options["software"] = option

    if reads:
        reads = check_path(reads)
        tgs_task, mito_reads = get_tgs_reads_task(
            prefix,
            database,
            reads,
            job_type,
            work_dir=os.path.join(work_dir, work_dict["tgs"]),
            out_dir=out_dir,
            thread=thread,
            gc=gc,
            base=base,
            score=score,
            sequencer=sequencer
        )
        dag.add_task(tgs_task)
    else:
        mito_reads = ""

    do_dag(dag, concurrent, refresh)

    return mito_r1, mito_r2, mito_reads, options


def get_mito(args):

    mito_r1, mito_r2, mito_reads, options = run_get_mito(
        prefix=args.prefix,
        database=args.database,
        platform=args.platform,
        read1=args.read1,
        read2=args.read2,
        sequencer=args.sequencer,
        reads=args.reads,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        thread=args.thread,
        gc=args.gc,
        base=args.base,
        score=args.score
    )

    with open(os.path.join(args.out_dir, "get_mito.json"), "w") as fh:
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
    get_mito Get mitochondrial chloroplast reads
version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_get_mito_args(parser)
    args = parser.parse_args()
    get_mito(args)


if __name__ == "__main__":
    main()
