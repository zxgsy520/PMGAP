#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
result structure


"""

from __future__ import absolute_import

import sys
import argparse
import logging
import json

from thirdparty.dagflow import Task, DAG, do_dag
from pmgap import __author__, __email__, __version__
from pmgap.config import *
from pmgap.common import mkdir, check_paths, get_version
from pmgap.parser import add_annotation_args
from pmgap.run_CDS import create_cds_annotation_dag


LOG = logging.getLogger(__name__)
__all__ = []


def create_prodigal_task(genome, prefix, work_dir=".", out_dir=".",
                         job_type="local", gcode=11, locus_tag="NPGAP",
                         meta=False):

    if meta:
        p = "meta"
    else:
        p = "single"

    option = OrderedDict()
    option["prodigal"] = {
        "version": get_version(SOFTWARE_VERSION["prodigal"]),
        "option": "-p %s -g %s" % (meta, gcode)
    }

    task = Task(
        id="prodigal",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={prodigal}:$PATH
prodigal -f gff -g {gcode} -p {p} -i {genome} -o {prefix}.prodigal.gff3
export PATH={python}:$PATH
python {script}/merge_gene.py {prefix}.prodigal.gff3 \\
  --locus {locus} --gcode {gcode} > {prefix}.merge.gff3
python {script}/get_gene.py --gff {prefix}.merge.gff3 --genome {genome} \\
  --type CDS --translate --gcode {gcode} > {prefix}.protein.fasta
""".format(prodigal=PRODIGAL_BIN,
            python=PYTHON_BIN,
            script=SCRIPTS,
            gcode=gcode,
            locus=locus_tag,
            p=p,
            out=out_dir,
            genome=genome,
            prefix=prefix)
    )
    protein = os.path.join(work_dir, "%s.protein.fasta" % prefix)
    gff3 = os.path.join(work_dir, "%s.merge.gff3" % prefix)

    return task, option, protein, gff3


def create_vgas_task(genome, prefix, work_dir=".", out_dir=".",
                         job_type="local", gcode=11, locus_tag=""):

    
    task = Task(
        id="vgas",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={vgas}:{python}:$PATH
rm -f blastp handle.py
ln -s {vgas}/blastp
ln -s {vgas}/handle.py
{vgas}/vgas {genome} {prefix} -i 1 -n -p -b
python {script}/process_vgas.py {prefix} --gene gene_{prefix} \\
  --protein protein_{prefix} --genome {genome} --prefix {prefix} \\
  --locus {locus} --gcode {gcode} > {prefix}.merge.gff3
mv {prefix}.protein.fasta {prefix}.genomic.pep
cut -f 1 {prefix}.genomic.pep >{prefix}.protein.fasta
""".format(vgas=VGAS_BIN,
            python=PYTHON_BIN,
            script=SCRIPTS,
            gcode=gcode,
            locus=locus_tag,
            out=out_dir,
            genome=genome,
            prefix=prefix)
    )
    protein = os.path.join(work_dir, "%s.protein.fasta" % prefix)
    gff3 = os.path.join(work_dir, "%s.merge.gff3" % prefix)

    return task, protein, gff3


def create_genemarks_task(genome, kingdom, prefix, work_dir=".", out_dir=".",
                         job_type="local", gcode=1, locus_tag=""):
    
    if kingdom == "viruses":
        x = "--virus"
        gcode = 1
    elif kingdom == "phage":
        x = "--phage"
        gcode = 11
    else:
        x = "--prok"
        gcode = 11

    task = Task(
        id="genemarks",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={genemarks}:{python}:$PATH
#cp {genemarks}/gm_key_64 ~/.gm_key
{genemarks}/gmsn.pl {x} --format GFF --fnn --faa --output {prefix}.mgm {genome}
python {script}/gmsn2gff.py {prefix}.mgm --locus {locus} --gcode {gcode} >{prefix}.merge.gff3
python {script}/get_gene.py --gff {prefix}.merge.gff3 --genome {genome} \\
  --type CDS --translate --gcode {gcode} > {prefix}.protein.fasta
""".format(genemarks=GENEMARKS_BIN,
            python=PYTHON_BIN,
            script=SCRIPTS,
            x=x,
            gcode=gcode,
            locus=locus_tag,
            out=out_dir,
            genome=genome,
            prefix=prefix)
    )
    protein = os.path.join(work_dir, "%s.protein.fasta" % prefix)
    gff3 = os.path.join(work_dir, "%s.merge.gff3" % prefix)

    return task, protein, gff3


def run_annotation(genome, prefix, locus_tag="NPGAP", organism="", strain="", template="",
                   kingdom="Bacteria", gcode=11, threads=1, evalue=1e-06, coverage=30,
                   job_type="local", refresh=30, concurrent=10, work_dir=".", out_dir="."):
    if not locus_tag:
        locus_tag = prefix
    LOG.info("ANNOTATE-1: gene predict\n")
    dag = DAG("gene_predict")

    genome = check_paths(genome)
    work_dict = {
        "stru": "Structure",
        "func": "Function",
    }
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        if "Structure" in v:
            continue
        mkdir(os.path.join(out_dir, v))

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    if kingdom != "viruses" and kingdom != "viral":
        prodigal_task, option, protein, gff3 = create_prodigal_task(
            genome=genome,
            prefix=prefix,
            gcode=gcode,
            locus_tag=locus_tag,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["stru"]),
            out_dir=out_dir
        )
        options["software"].update(option)
        dag.add_task(prodigal_task)
    else:
        #vgas_task, protein, gff3 = create_vgas_task(
        #    genome=genome,
        #    prefix=prefix,
        #    gcode=gcode,
        #    locus_tag=locus_tag,
        #    job_type=job_type,
        #    work_dir=os.path.join(work_dir, work_dict["stru"]),
        #    out_dir=out_dir
        #)
        #dag.add_task(vgas_task)
        genemarks_task, protein, gff3 = create_genemarks_task(
            genome=genome,
            kingdom=kingdom,
            prefix=prefix,
            work_dir=os.path.join(work_dir, work_dict["stru"]),
            out_dir=out_dir,
            job_type=job_type,
            gcode=gcode,
            locus_tag=locus_tag
        )
        dag.add_task(genemarks_task)

    do_dag(dag, refresh_time=refresh, concurrent_tasks=concurrent)

    LOG.info("ANNOTATE-2: gene annotate\n")
    dag = DAG("gene_annotate")

    cds_dag, option = create_cds_annotation_dag(
        protein=protein,
        prefix=prefix,
        kingdom=kingdom,
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["func"]),
        out_dir=os.path.join(out_dir, work_dict["func"])
    )

    options["software"].update(option["software"])
    options["database"].update(option["database"])

    tmp = ""
    if organism:
        tmp += " --organism '%s'" % organism
    if strain:
        tmp += " --strain '%s'" % strain
    tmp2 = ""
    if template:
        tmp2 = "-t %s" % check_paths(template)

    options["software"]["tbl2asn"] = {
        "version": get_version(SOFTWARE_VERSION["tbl2asn"]),
        "option": "default"
    }

    merge = Task(
        id="merge_annotation",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={tbl2asn}:{python}:$PATH
python {script}/update_annotation.py {gff3} --refseq {cds}/{prefix}.refseq.tsv \\
  --swissprot {cds}/{prefix}.SwissProt.tsv --COG {cds}/{prefix}.COG.tsv \\
  --TIGRFAMs {cds}/{prefix}.TIGRFAMs.tsv --Pfam {cds}/{prefix}.Pfam.tsv \\
  --KEGG {cds}/{prefix}.KEGG.tsv --GO {cds}/{prefix}.GO.tsv -p {prefix} >{prefix}.genomic.gff3

python {script}/merge_stat_anno.py {gff3} --refseq {cds}/{prefix}.refseq.tsv \\
  --swissprot {cds}/{prefix}.SwissProt.tsv --COG {cds}/{prefix}.COG.tsv \\
  --TIGRFAMs {cds}/{prefix}.TIGRFAMs.tsv --Pfam {cds}/{prefix}.Pfam.tsv \\
  --KEGG {cds}/{prefix}.KEGG.tsv --GO {cds}/{prefix}.GO.tsv >{prefix}.merge.annotate.xls

python {script}/gff2tbl_bac.py {prefix}.genomic.gff3 > {prefix}.genomic.tbl
python {script}/process_assembly.py {genome} --gcode {gcode} {args} > {prefix}.genomic.fasta
tbl2asn -i {prefix}.genomic.fasta -V b -s T {args2}

mv {prefix}.genomic.gbf {prefix}.genomic.gb
cp {prefix}.merge.annotate.xls {prefix}.genomic.gff3 {prefix}.genomic.gb {prefix}.genomic.sqn {prefix}.function_summary.tsv {out}
cd {out}

python {script}/get_gene.py --gff {prefix}.genomic.gff3 --genome {genome} --type CDS --translate --gcode {gcode} > {prefix}.protein.fasta
python {script}/get_gene.py --gff {prefix}.genomic.gff3 --genome {genome} --type CDS > {prefix}.RNA.fasta
python {script}/cds_stat.py {prefix}.RNA.fasta --name {prefix} --min_length 0 --gcode {gcode}
python {script}/gene_stat.py --gff {prefix}.genomic.gff3 --fasta {genome} > {prefix}.structure_summary.tsv

""".format(tbl2asn=TBL2ASN_BIN,
            python=PYTHON_BIN,
            script=SCRIPTS,
            prefix=prefix,
            genome=genome,
            gff3=gff3,
            gcode=gcode,
            args=str(tmp),
            args2=tmp2,
            out=out_dir,
            cds=os.path.join(out_dir, work_dict["func"]),
           )
    )

    merge.set_upstream(*cds_dag.tasks.values())
    dag.add_task(*cds_dag.tasks.values())
    dag.add_task(merge)

    with open(os.path.join(out_dir, "annotate.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    do_dag(dag, refresh_time=refresh, concurrent_tasks=concurrent)

    return os.path.join(out_dir, "%s.genomic.gff3" % prefix), os.path.join(out_dir, "%s.protein.fasta" % prefix), os.path.join(out_dir, "%s.RNA.fasta" % prefix)


def annotate(args):

    run_annotation(
        genome=args.genome,
        prefix=args.prefix,
        locus_tag=args.locustag,
        kingdom=args.kingdom,
        organism=args.organism,
        strain=args.strain,
        template=args.template,
        gcode=args.gcode,
        threads=args.threads,
        evalue=args.evalue,
        coverage=args.coverage,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir)


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
NPGAP CDS annotation pipeline

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_annotation_args(parser)

    annotate(parser.parse_args())


if __name__ == "__main__":
    main()
