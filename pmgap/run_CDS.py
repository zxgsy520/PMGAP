#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
*.refseq.m6      blast results of refseq database
*.refseq.tsv     refseq annotation file
*.KEGG.m6        blast results of KEGG database
*.KEGG.tsv       KEGG annotation file
*.COG.out        blast results of COG database
*.COG.tsv        COG annotation file
*.ipr.out        interproscan annotation file
*.ipr.tsv        interproscan annotation file
"""

import sys
import argparse
import logging

from thirdparty.dagflow import Task, ParallelTask, DAG, do_dag
from pmgap import __author__, __email__, __version__
from pmgap.config import *
from pmgap.common import check_paths, cd, mkdir, get_version
from pmgap.parser import add_cds_args
from thirdparty.seqkit.split import seq_split


LOG = logging.getLogger(__name__)
__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def create_refseq_task(proteins, prefix, db, evalue, coverage, threads,
                       job_type, work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "Refseq"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option= "-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue 1e-05 --threads {threads} --out {{prefixs}}.refseq.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           script=SCRIPTS,
           threads=threads),
        prefixs=prefixs,
        protein=proteins,
    )

    join_task = Task(
        id="merge_refseq",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.refseq.m6 > {prefix}.refseq.m6
time {script}/blast_filter.py {prefix}.refseq.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend stitle evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.refseq.out
{script}/refseqproc.py {prefix}.refseq.out > {prefix}.refseq.tsv
{script}/sure_species.py -i {prefix}.refseq.tsv -n {prefix}

cp {prefix}.refseq.tsv {prefix}.refseq.m6 {prefix}.species.tsv {out_dir}
""".format(id=id,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_interproscan_task(proteins, prefix, job_type, work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id="ipr"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={interproscan}:$PATH
time interproscan.sh -i {{protein}} -appl Pfam,TIGRFAM,SMART -iprlookup -goterms -t p -f TSV -o {{prefixs}}.ipr.out
""".format(interproscan=INTERPROSCAN_BIN,
           script=SCRIPTS,
           ),
        protein=proteins,
        prefixs=prefixs

    )

    join_task = Task(
        id="merge_ipr",
        work_dir=work_dir,
        type="local",
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={script}:$PATH
cat {id}*/*.ipr.out > {prefix}.ipr.out
time iprproc.py {prefix}.ipr.out --tigrfams {tigrfams}/TIGRFAMS.link --prefix {prefix}  > {prefix}.other.tsv

goproc.py --go {prefix}.WEGO.txt --go_obo {go}/go.obo > {prefix}.GO.tsv
new.get_GO.classify.py {prefix}.WEGO.txt {go}/GO_level4.deal.txt > {prefix}.go.classify.xls
{rscript} {script}/get.level2_3.counts.R {prefix}.go.classify.xls {prefix}.go2.xls {prefix}.go3.xls
#{rscript} {script}/GO.level.stat.R {prefix}.go2.xls {prefix}.WEGO.pdf
#ghostscript -dNOSAFER -r600 -dBATCH -sDEVICE=pngalpha -dNOPAUSE -dEPSCrop -sOutputFile={prefix}.WEGO.png {prefix}.WEGO.pdf
#cp *.pdf *.png {out_dir}
cp {prefix}.ipr.out {prefix}.*.tsv {out_dir}
""".format(id=id,
           prefix=prefix,
           tigrfams=TIGRFAMS,
           go=GO,
           rscript=RSCRIPT,
           script=SCRIPTS,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_kegg_task(proteins, prefix, db, evalue, coverage, threads, job_type,
                     work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "KEGG"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue {evalue} --threads {threads} --out {{prefixs}}.KEGG.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           threads=threads),
        protein=proteins,
        prefixs=prefixs,
    )

    join_task = Task(
        id="merge_KEGG",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.KEGG.m6 > {prefix}.KEGG.m6
time {script}/blast_filter.py {prefix}.KEGG.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.KEGG.out

{script}/keggproc.py {prefix}.KEGG.out --ko {pep2ko} --pathway {pathway} --keg {keg} > {prefix}.KEGG.tsv
{script}/get_kegg_pathway.py --input {prefix}.KEGG.out -pd {pep2ko} \\
--pathway {pathway} --keg {keg} --out {prefix}.kegg_pathway.tsv
cut -f 1,3,4 {prefix}.KEGG.tsv > {prefix}.KEGG.KO
{script}/make_keg.py --keg {keg} --in {prefix}.KEGG.KO --out {prefix}.KEGG --plot

#cp {prefix}.KEGG.png {prefix}.KEGG.pdf {out_dir}
cp {prefix}.kegg_pathway.tsv {prefix}.KEGG.KO {prefix}.KEGG.keg {prefix}.KEGG.tsv {prefix}.KEGG.m6 {out_dir}
""".format(id=id,
           pep2ko=KEGG_KO,
           pathway=KEGG_PATH,
           keg=KEGG_KEG,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_cog_task(proteins, prefix, db, evalue, coverage, threads, job_type,
                    work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "COG"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue {evalue} --threads {threads} --out {{prefixs}}.COG.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           threads=threads),
        protein=proteins,
        prefixs=prefixs
    )

    join_task = Task(
        id="merge_COG",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.COG.m6 > {prefix}.COG.m6
time {script}/blast_filter.py {prefix}.COG.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend slen stitle evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.COG.out
time {script}/cogproc.py {prefix}.COG.out --name {cog} > {prefix}.COG.tsv
{script}/plot_cog.py {prefix}.COG.tsv --func {func} --name {cog} --out {prefix}.COG
#cp {prefix}.COG.pdf {prefix}.COG.png {out_dir}
cp {prefix}.COG.m6 {prefix}.COG.tsv {out_dir}
""".format(id=id,
           prefix=prefix,
           cog=COG_NAME,
           func=COG_FUNC,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_swissprot_task(proteins, prefix, db, evalue, coverage, threads, job_type,
                    work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "SwissProt"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue {evalue} --threads {threads} --out {{prefixs}}.SwissProt.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           threads=threads),
        protein=proteins,
        prefixs=prefixs
    )

    join_task = Task(
        id="merge_SwissProt",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.SwissProt.m6 > {prefix}.SwissProt.m6
time {script}/blast_filter.py {prefix}.SwissProt.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend slen stitle evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.SwissProt.out
time {script}/swissprotproc.py {prefix}.SwissProt.out > {prefix}.SwissProt.tsv
cp {prefix}.SwissProt.m6 {prefix}.SwissProt.tsv {out_dir}
""".format(id=id,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_cds_annotation_dag(protein, prefix, kingdom, evalue, coverage,
                              threads, job_type, work_dir, out_dir):

    kingdom = kingdom.lower()
    dag = DAG("CDS")

    work_dict = {
        "split": "00_data",
        "refseq": "01_Refseq",
        "ipr": "02_ipr",
        "kegg": "03_KEGG",
        "cog": "04_COG",
        "swissprot": "05_SwissProt"
    }

    work_dir = mkdir(os.path.abspath(work_dir))
    out_dir = mkdir(os.path.abspath(out_dir))
    protein = check_paths(protein)

    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))

    proteins = seq_split([protein], mode="length", num=5000000, output_dir=work_dict["split"])

    _options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    _options["software"]["blastp"] = {
        "version": get_version(SOFTWARE_VERSION["blastp"]),
        "option": "-evalue %s -outfmt '6 std qlen slen stitle' -max_target_seqs 5" % evalue
    }
    _options["database"]["refseq"] = {
        "version": os.path.basename(REFSEQ),
    }
    refseq_tasks, refseq_join = create_refseq_task(
        proteins=proteins,
        prefix=prefix,
        db=REFSEQ_TAXON[kingdom],
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["refseq"],
        out_dir=out_dir
    )
    dag.add_task(*refseq_tasks)
    dag.add_task(refseq_join)

    _options["software"]["Interproscan"] = {
        "version": get_version(SOFTWARE_VERSION["interproscan"]),
        "option": "-appl Pfam,TIGRFAM,SMART -iprlookup -goterms -t p -f TSV"
    }
    _options["database"]["Pfam"] = {
        "version": os.path.basename(PFAM),
    }
    _options["database"]["TIGRFAMs"] = {
        "version": os.path.basename(TIGRFAMS),
    }
    _options["database"]["GO"] = {
        "version": os.path.basename(GO),
    }
    ipr_tasks, ipr_join = create_interproscan_task(
        proteins=proteins,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dict["ipr"],
        out_dir=out_dir
    )
    dag.add_task(*ipr_tasks)
    dag.add_task(ipr_join)


    _options["database"]["KEGG"] = {
        "version": os.path.basename(KEGG),
    }
    kegg_tasks, kegg_join = create_kegg_task(
        proteins=proteins,
        prefix=prefix,
        db=KEGG_TAXON[kingdom],
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["kegg"],
        out_dir=out_dir
    )
    dag.add_task(*kegg_tasks)
    dag.add_task(kegg_join)


    _options["software"]["rpsblast"] = {
        "version": get_version(SOFTWARE_VERSION["rpsblast"]),
        "option": "-evalue 0.01 -seg no -outfmt 5"
    }
    _options["database"]["COG"] = {
        "version": os.path.basename(COG),
    }
    cog_tasks, cog_join = create_cog_task(
        proteins=proteins,
        prefix=prefix,
        db=COG_TAXON[kingdom],
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["cog"],
        out_dir=out_dir
    )
    dag.add_task(*cog_tasks)
    dag.add_task(cog_join)

    swissprot_tasks, swissprot_join = create_swissprot_task(
        proteins=proteins,
        prefix=prefix,
        db=SWISSPROT_DB,
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["swissprot"],
        out_dir=out_dir
    )
    dag.add_task(*swissprot_tasks)
    dag.add_task(swissprot_join)

    return dag, _options


def run_cds_annotation(proteins, prefix, kingdom, evalue, coverage, threads, job_type, concurrent, refresh, work_dir, out_dir):

    dag, options = create_cds_annotation_dag(proteins, prefix, kingdom, evalue, coverage, threads, job_type, work_dir, out_dir)
    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
CDS annotation pipeline

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_cds_args(parser)
    args = parser.parse_args()

    run_cds_annotation(
        proteins=args.protein,
        prefix=args.name,
        kingdom=args.kingdom,
        evalue=args.evalue,
        coverage=args.coverage,
        threads=args.threads,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )


if __name__ == "__main__":
    main()
