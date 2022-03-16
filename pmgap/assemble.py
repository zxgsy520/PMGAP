#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import argparse
import logging

from thirdparty.dagflow import DAG, Task, do_dag
from pmgap.config import *
from pmgap.parser import add_assemble_args
from pmgap.common import check_path, mkdir, get_version

LOG = logging.getLogger(__name__)

__version__ = "v1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def create_unicycler_task(prefix, read1, read2, reads, job_type, work_dir, out_dir,
                         minlen=500, thread=4):

    temp = ""
    if reads:
        reads = "--long %s" % reads
        temp = "--long %s"
    option = {}
    option["unicycler"] = {
        "version": get_version(SOFTWARE_VERSION["unicycler"]),
        "option": "-1 read1 -2 read2 %s --keep 3 --mode normal" % temp
    }
    option["bwa"] = {
        "version": get_version(SOFTWARE_VERSION["bwa"]),
        "option:": "default"
    }
    option["circlator"] = {
        "version": get_version(SOFTWARE_VERSION["circlator"]),
        "option": "circlator fixstart",
    }

    task=Task(
        id='unicycler',
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={unicycler}:$PATH
unicycler -1 {read1} -2 {read2} \\
  {reads} --keep 3 --mode normal -o {prefix} -t {thread}
python {script}/add_circular_edge.py {prefix}/assembly.fasta --minlen {minlen} >{prefix}.contigs.fasta
mv {prefix}/assembly.gfa {prefix}.assembly.gfa
#rm -rf {prefix}
""".format(unicycler=UNICYCLER_BIN,
            script=SCRIPTS,
            prefix=prefix,
            read1=read1,
            read2=read2,
            reads=reads,
            minlen=minlen,
            thread=thread,
        )
    )

    circlator_task = Task(
        id='circular',
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={python}:$PATH
python {script}/check_circularity.py --blastn {blast}/blastn {prefix}.contigs.fasta > {prefix}.circularise.fasta
python {script}/circlatorproc.py {prefix}.circularise.fasta TEMP

if [ -e TEMP.circular.fasta ]; then
    export PATH={prodigal}:{mummer}:{circlator}:$PATH
    circlator fixstart TEMP.circular.fasta TEMP
    cat TEMP.linear.fasta TEMP.fasta > {prefix}.nr.fasta
else
    cat TEMP.linear.fasta > {prefix}.nr.fasta
fi
export PATH={python}:$PATH
python {script}/process_assembly.py {prefix}.nr.fasta --format "tig%05d" --min_length 0 > {prefix}.genomic.fasta
python {script}/genome_stat.py {prefix}.genomic.fasta --contig {prefix}.contig_stat.tsv --base {prefix}.base_stat.tsv

cp {prefix}.genomic.fasta {prefix}.contig_stat.tsv {prefix}.base_stat.tsv {out_dir}
rm TEMP.*
""".format(python=PYTHON_BIN,
            blast=BLAST_BIN,
            prodigal=PRODIGAL_BIN,
            mummer=MUMMER_BIN,
            circlator=CIRCLATOR_BIN,
            script=SCRIPTS,
            prefix=prefix,
            out_dir=out_dir
        )
    )
    task.set_downstream(circlator_task)

    return task, circlator_task, option, os.path.join(out_dir, '%s.genomic.fasta' % prefix)


def create_gc_depth_task(genome, prefix, read1, read2, job_type,
                         work_dir, out_dir, window=100, thread=4):

    option = {}
    option["minimap2"] = {
        "version": get_version(SOFTWARE_VERSION["minimap2"]),
        "option": "default"
    }
    option["samtools"] = {
        "version": get_version(SOFTWARE_VERSION["samtools"]),
        "option:": "default"
    }

    task = Task(
        id="depth",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={samtools}:{minimap2}:{python}:$PATH
if [ ! -e bam_done ]; then
    minimap2 -t {thread} -ax sr {genome} {read1} {read2} |samtools view -bS |samtools sort -o {prefix}.sort.bam
    samtools depth -aa {prefix}.sort.bam > {prefix}.depth
    touch bam_done
    rm {prefix}.sort.bam
fi
python {script}/stat_length_gc.py {genome} -d {prefix}.depth >{prefix}.length_gc.xls
python {script}/plot_depth_stat.py {prefix}.depth --window {window} --prefix {prefix}
cp {prefix}.length_gc.xls {prefix}.depth.png {prefix}.depth.pdf {out_dir}
""".format(samtools=SAMTOOLS_BIN,
            minimap2=MINIMAP2_BIN,
            script=SCRIPTS,
            python=PYTHON_BIN,
            genome=genome,
            read1=read1,
            read2=read2,
            prefix=prefix,
            thread=thread,
            window=window,
            out_dir=out_dir
        )
    )

    return task, option, os.path.join(work_dir, "%s.gc_depth.png" % prefix)


def run_assembly(prefix, read1, read2, reads, thread, job_type,
                 concurrent, refresh, work_dir, out_dir):

    read1 = check_path(read1)
    read2 = check_path(read2)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    dag = DAG("assemble")
    unicycler_task, circlator_task, option, genome = create_unicycler_task(
        prefix=prefix,
        read1=read1,
        read2=read2,
        reads=reads,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir,
        thread=thread,
        minlen=500
    )
    options["software"].update(option)

    gc_depth_task, option, gc_depth = create_gc_depth_task(
        genome=genome,
        read1=read1,
        read2=read2,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir,
        thread=thread,
        window=100
    )
    options["software"].update(option)

    circlator_task.set_downstream(gc_depth_task)

    dag.add_task(unicycler_task)
    dag.add_task(circlator_task)
    dag.add_task(gc_depth_task)

    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    return genome, options


def assemble(args):

    genome, options = run_assembly(
        prefix=args.prefix,
        read1=args.read1,
        read2=args.read2,
        thread=args.thread,
        job_type=args.job_type,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        concurrent=args.concurrent,
        refresh=args.refresh
    )
    with open(os.path.join(args.out_dir, "assemble.json"), "w") as fh:
        json.dump(options, fh, indent=2)


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

    parser = add_assemble_args(parser)
    args = parser.parse_args()
    assemble(args)


if __name__ == "__main__":
    main()
