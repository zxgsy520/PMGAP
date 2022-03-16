#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import os.path
import argparse
import logging
import shutil

from pmgap.config import *
from pmgap.common import mkdir, check_path, check_paths, rm #read_config
from pmgap.ngs_qc import run_ngs_qc
from pmgap.assemble import run_assembly
from pmgap.get_mito import run_get_mito
from pmgap.annotate import run_annotation
from pmgap.abricate import run_abricate
from pmgap.backup import backup
from pmgap.parser import add_all_args

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def run_report(prefix, project, id, platform, author, reviewer, out_dir, work_dir):

    shell = """
{script}/report_ngs.py {work_dir}/report_config.cfg \\
  --fasta {result}/02_Assembly/{name}.genomic.fasta \\
  --reads {result}/01_Data/{name}.qc.xls \\
  --contig {result}/02_Assembly/{name}.contig_stat.tsv \\
  --base {result}/02_Assembly/{name}.base_stat.tsv \\
  --length_gc {result}/02_Assembly/{name}.length_gc.xls \\
  --codon {result}/03_Annotation/{name}.codon.tsv \\
  --structure {result}/03_Annotation/{name}.structure_summary.tsv \\
  --function {result}/03_Annotation/{name}.function_summary.tsv \\
  --jsons {result}/pmgap.json  {result}/03_Annotation/annotate.json \\
  --base_content {result}/01_Data/{name}.base_content.png \\
  --base_gc {result}/01_Data/{name}.base_gc.png \\
  --base_quality {result}/01_Data/{name}.base_quality.png \\
  --depth {result}/02_Assembly/{name}.depth.png \\
  --protein {result}/03_Annotation/{name}.protein_length.png \\
  --html {template}/ngs_html \\
  --out {out_dir}
""".format(script=SCRIPTS,
        name=prefix,
        template=TEMPLATES,
        result=os.path.join(out_dir, "result"),
        work_dir= work_dir,
        out_dir=out_dir,
    )

    config ="""[general]
project={project}
id={id}
name={name}
sequencer={platform}
author={author}
reviewer={reviewer}
lib_num=1
cell_num=1
""".format(name=prefix,
        project=project,
        id=id,
        platform=platform,
        author=author,
        reviewer=reviewer)

    run = os.path.join(work_dir, "report.sh")
    fs = open(run, 'w')
    fc = open(os.path.join(work_dir, "report_config.cfg"), 'w')
    fs.write(shell)
    fc.write(config)
    fs.close()
    fc.close()

    os.system("sh %s" % run)


def run_all(prefix, taxon, platform, reads1, reads2, sequencer, reads,
            database, project, projectid, job_type, work_dir, out_dir,
            concurrent=10, refresh=15, trim=3, thread=4, gc=20,
            base="all", score=0, gcode=11):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    reads1 = check_paths(reads1)
    reads2 = check_paths(reads2) 

    work_dict = {
        "data": "01_Data",
        "asm": "02_Assembly",
        "ann": "03_Annotation",
    }
    
    result_dir = mkdir(os.path.join(out_dir, "result"))

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(result_dir, v))

    LOG.info("run qc")
    clean1, clean2, stat_qc, quality, content, base_gc, options = run_ngs_qc(
        reads1=reads1,
        reads2=reads2,
        prefix=prefix,
        trim=trim,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=os.path.join(work_dir, work_dict["data"]),
        out_dir=os.path.join(result_dir, work_dict["data"])
    )
    if taxon in ["mitochondrion", "plastid", "viruses"]:
        if taxon in ["viruses"]:
           gc = 0
        LOG.info("run get_mito")
        if database:
            database = check_path(database)
            LOG.info("database:%s" %  database)
        #else:
        #    database = check_path(REFSEQ_TAXON[taxon])
        #    LOG.info("database:%s" %  database)
            clean1, clean2, reads, options = run_get_mito(
                prefix=prefix,
                database=database,
                platform=platform,
                read1=clean1,
                read2=clean2,
                sequencer=sequencer,
                reads=reads,
                job_type=job_type,
                concurrent=concurrent,
                refresh=refresh,
                work_dir=os.path.join(work_dir, work_dict["asm"]),
                out_dir=os.path.join(result_dir, work_dict["asm"]),
                thread=thread,
                gc=gc,
                base=base,
                score=score
            )
        else:
            pass
     
    LOG.info("run assembly") 
    genome, options = run_assembly(
        prefix=prefix,
        read1=clean1,
        read2=clean2,
        reads=reads,
        thread=thread,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["asm"]),
        out_dir=os.path.join(result_dir, work_dict["asm"]),
        concurrent=concurrent,
        refresh=refresh
    )
    
    LOG.info("run Annotation")
    gff, protein, gene = run_annotation(
        genome=genome,
        prefix=prefix,
        locus_tag="",
        kingdom=taxon,
        organism="Unknow",
        strain=prefix,
        template="",
        gcode=gcode,
        threads=thread,
        evalue=1e-05,
        coverage=30,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=os.path.join(work_dir, work_dict["ann"]),
        out_dir=os.path.join(result_dir, work_dict["ann"])
    )

    if "plasmid" in taxon:
        LOG.info("run abricate")
        run_abricate(gene=gene,
            prefix=prefix,
            threads=thread,
            job_type=job_type,
            concurrent=concurrent,
            refresh=refresh,
            work_dir=mkdir(os.path.join(work_dir, "04_Abricate")),
            out_dir=mkdir(os.path.join(result_dir, "04_Abricate")),
            database="card;vfdb"
        )

    with open(os.path.join(result_dir, "pmgap.json"), "w") as fh:
        json.dump(options, fh, indent=2)
    
    run_report(
        prefix=prefix,
        project=project,
        id=projectid,
        platform=platform,
        author="百易汇能",
        reviewer="百易汇能",
        out_dir=out_dir,
        work_dir=work_dir)

    LOG.info("run_backup")
    backup(
        name=prefix,
        project=projectid,
        job_type="local",
        work_dir=work_dir,
        out_dir=out_dir
    )

    return 0


def _all(args):

    run_all(prefix=args.prefix,
            taxon=args.taxon,
            platform=args.platform,
            reads1=args.reads1,
            reads2=args.reads2,
            sequencer=args.sequencer,
            reads=args.reads,
            database=args.database,
            project=args.project,
            projectid=args.projectid,
            job_type=args.job_type,
            work_dir=args.work_dir,
            out_dir=args.out_dir,
            concurrent=args.concurrent,
            refresh=args.refresh, 
            trim=args.trim,
            thread=args.thread,
            gc=args.gc,
            base=args.base,
            score=args.score,
            gcode=args.gcode,
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

    parser = add_all_args(parser)
    args = parser.parse_args()
    
    _all(args)


if __name__ == "__main__":
    main()
