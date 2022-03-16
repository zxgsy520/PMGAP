#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../"))
from pmgap.config import *
from pmgap.common import check_path, check_paths, mkdir, read_tsv
from thirdparty.dagflow import DAG, Task, do_dag

LOG = logging.getLogger(__name__)
__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []

try:
    reload(sys)
    sys.setdefaultencoding('utf-8')
except:
    pass

def create_pmgap_task(prefix, reads1, reads2, project, projectid,
                      job_type, work_dir, out_dir,
                      taxon, platform="mgi", database="", gcode=11):

    if database:
       database = check_path(database)
       database = "--database %s" % database
    task = Task(
        id="pmgap_%s" % prefix,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{root}/pmgap.py all \\
  --reads1 {reads1} \\
  --reads2 {reads2} \\
  --prefix {prefix} --taxon {taxon} --trim 5 --thread 10 --job_type {job_type} \\
  --platform {platform} {database} \\
  --project {project} --projectid {projectid} \\
  --work_dir {work}/{prefix}  --out_dir {out}/{prefix} --gcode {gcode}
""".format(root=ROOT,
            prefix=prefix,
            reads1=" ".join(reads1),
            reads2=" ".join(reads2),
            taxon=taxon,
            gcode=gcode,
            platform=platform,
            database=database,
            project=project,
            projectid=projectid,
            job_type=job_type,
            work=work_dir,
            out=out_dir
        )
    )

    return task


def run_pmgap_multi(data, taxon, project, projectid,
                    work_dir, out_dir, concurrent, refresh,
                    job_type="local", platform="mgi", database="", gcode=11):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    data = check_path(data)
    r = {}

    for line in read_tsv(data, "\t"):
        if line[0] not in r:
             r[line[0]] = [[], []]
        r[line[0]][0].append(line[1])
        r[line[0]][1].append(line[2])

    dag = DAG("run_pmgap_multi")
    for i in r:
        reads1 = check_paths(r[i][0])
        reads2 = check_paths(r[i][1])
        task = create_pmgap_task(
            prefix=i,
            reads1=reads1,
            reads2=reads2,
            taxon=taxon,
            gcode=gcode,
            platform=platform,
            database=database,
            project=project,
            projectid=projectid,
            job_type=job_type,
            work_dir=work_dir,
            out_dir=out_dir
        )
        LOG.info("run %s" % i)
        dag.add_task(task)

    do_dag(dag, concurrent, refresh)

    return 0


def pmgap_multi(args):

    run_pmgap_multi(
        data=args.data,
        taxon=args.taxon,
        gcode=args.gcode,
        project=args.project,
        projectid=args.projectid,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        concurrent=args.concurrent,
        refresh=args.refresh,
        job_type=args.job_type,
        platform=args.platform,
        database=args.database,
    )


def add_pmgap_multi_args(parser):

    parser.add_argument("data", metavar='FILE', type=str,
        help="Input sequencing list file.")
    parser.add_argument("--taxon", choices=["mitochondrion", "plasmid", "plastid", "viruses"], default="mitochondrion",
        help="Choose taxon for assembly species,default=plasmid.")
    parser.add_argument("-pl", "--platform", choices=["illumina", "mgi"], default="mgi",
        help="Select sequencing platform, default=mgi")
    parser.add_argument("-d", "--database", metavar="FILE", type=str, default="",
        help="Input the reference database path, default="".")
    parser.add_argument("--gcode", choices=GENETIC_CODES, type=int,
                        help="Genetic codes used, see "
                             "'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c'"
                             "for more information (default: 11)", default=11),
    parser.add_argument("-pro", "--project", metavar="STR", type=str, required=True,
        help="Input project name, default=None.")
    parser.add_argument("-pid", "--projectid", metavar="STR", type=str, required=True,
        help="Input project id, default=None.")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10).")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30).")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local).")
    parser.add_argument("--work_dir", metavar="DIR", type=str, default=".",
        help="Work directory (default: current directory).")
    parser.add_argument("--out_dir", metavar="DIR", type=str, default=".",
        help="Output directory (default: current directory).")

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
attention:
    pmgap_multi.py genomes.list
File format：
name    r1	r2

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_pmgap_multi_args(parser)
    args = parser.parse_args()
    pmgap_multi(args)


if __name__ == "__main__":
    main()
