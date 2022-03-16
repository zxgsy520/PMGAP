#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pmgap.config import GENETIC_CODES

__all__ = ["add_ngs_qc_args", "add_get_mito_args", "add_assemble_args",
           "add_annotation_args", "add_all_args",
           ]


def add_workflow_args(parser):
    """
    add workflow arguments to parser
    :param parser: argparse object
    :return: parser
    """

    workflow_group = parser.add_argument_group(title="Workflow arguments", )
    workflow_group.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    workflow_group.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    workflow_group.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    workflow_group.add_argument("--work_dir", metavar="DIR", default="NPGAP.work",
        help="Work directory (default: current directory)")
    workflow_group.add_argument("--out_dir", metavar="DIR", default="NPGAP.out",
        help="Output directory (default: current directory)")

    return parser



def add_ngs_qc_args(parser):

    parser.add_argument("-r1", "--reads1", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--reads2", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r2.")
    parser.add_argument("-p", "--prefix", metavar="FILE", required=True,
        help="Sample name.")
    parser.add_argument("--trim", metavar="INT", type=int, default=5,
        help="Set trim length, default=5")
    parser.add_argument("-t", "--thread", type=int, default=4,
        help="Set the number of threads to run.")
    parser = add_workflow_args(parser)

    return parser


def add_get_mito_args(parser):
    
    parser.add_argument("-p", "--prefix", metavar="FILE", required=True,
        help="Sample name.")
    parser.add_argument("-d", "--database", metavar="FILE", required=True,
        help="Input the reference database path.")
    parser.add_argument("-pl", "--platform", choices=["illumina", "mgi"], default="illumina",
        help="Select sequencing platform, default=illumina")
    parser.add_argument("-r1", "--read1", metavar="FILE", type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", type=str, required=True,
        help="Second generation sequencing of reads r2.")
    parser.add_argument("-s", "--sequencer", choices=["PromethION", "GridION", "RSII", "Sequel"], default="PromethION",
        help="Used sequencer, default=PromethION")
    parser.add_argument("-r", "--reads", metavar="FLIE", type=str, default="",
        help="Input the third-generation sequencing reads.")
    parser.add_argument("-t", "--thread", type=int, default=4,
        help="Set the number of threads to run.")
    parser.add_argument("--gc", metavar="INT", type=int, default=0,
        help="Set the maximum GC content retained, default=0")
    parser.add_argument("--base", metavar="STR", type=str, default="all",
        help="Set the number of reserved bases, default=all")
    parser.add_argument("--score", metavar="INT", type=int, default=0,
        help="Match score, default=0")

    parser = add_workflow_args(parser)

    return parser


def add_assemble_args(parser):
   
    parser.add_argument("-p", "--prefix", metavar="FILE", required=True,
        help="Sample name.")
   # parser.add_argument("--taxon", choices=["mitochondrion", "plasmid", "plastid", "viral"], default="mitochondrion",
    #    help="Choose taxon for assembly species,default=plasmid.")
   # parser.add_argument("-pl", "--platform", choices=["illumina", "mgi"], default="illumina",
   #     help="Select sequencing platform, default=illumina")
    parser.add_argument("-r1", "--read1", metavar="FILE", type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", type=str, required=True,
        help="Second generation sequencing of reads r2.")
   # parser.add_argument("-s", "--sequencer", choices=["PromethION", "GridION", "RSII", "Sequel"], default="PromethION",
   #     help="Used sequencer, default=PromethION")
    parser.add_argument("-r", "--reads", metavar="FLIE", type=str, default="",
        help="Input the third-generation sequencing reads.")
#    parser.add_argument("-d", "--database", metavar="FILE", default="",
   #     help="Input the reference database path.")
 #   parser.add_argument("--gc", metavar="INT", type=int, default=0,
  #      help="Set the maximum GC content retained, default=0")
  #  parser.add_argument("--base", metavar="STR", type=str, default="all",
   #     help="Set the number of reserved bases, default=all")
   # parser.add_argument("--score", metavar="INT", type=int, default=0,
    #    help="Match score, default=0")
    parser.add_argument("-t", "--thread", type=int, default=4,
        help="Set the number of threads to run.")
    parser = add_workflow_args(parser)

    return parser


def add_cds_args(parser):
    """
    CDS annotation arguments
    :param parser:
    :return: parser
    """

    parser.add_argument("-p", "--protein", required=True, help="Proteins in fasta file")
    parser.add_argument("-n", "--name", metavar="STR", required=True,
                        help="Name used in pipeline")
    parser.add_argument("-l", "--locustag", metavar="STR", required=True,
                        help="Prefix of gene id")

    parser.add_argument("--kingdom", choices=["Archaea", "Bacteria", "Mitochondria", "Viruses"],
                        help="Kingdom of the sample (default: Bacteria)", default="Bacteria"),
    parser.add_argument("--gcode", choices=GENETIC_CODES, type=int,
                        help="Genetic codes used, see "
                             "'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c'"
                             "for more information (default: 11)", default=11),
    parser.add_argument("--meta", action="store_true",
                        help="True if sample is meta or plasmid")

    parser = add_blast_args(parser)
    parser = add_workflow_args(parser)

    return parser


def add_blast_args(parser):
    """
    add blast arguments to parser
    :param parser: argparse object
    :return: parser
    """
    parser.add_argument("--threads", metavar="INT", type=int,
                        help="Threads used to run blastp (default: 1)", default=1)
    parser.add_argument("--evalue", metavar="NUM", type=float,
                        help="Evalue cutoff of blastp for Refseq and KEGG (default: 1e-05)", default=1e-05, )
    parser.add_argument("--coverage", metavar="NUM", type=float,
                        help="Coverage cutoff of blastp for Refseq and KEGG (default: 30)", default=30,)

    return parser


def add_annotation_args(parser):

    parser.add_argument("-g", "--genome", metavar="FILE", required=True,
                        help="Genome in fasta file")
    parser.add_argument("-p", "--prefix", metavar="STR", required=True,
                        help="Name used in pipeline")
    parser.add_argument("-l", "--locustag", metavar="STR", required=True,
                        help="Prefix of gene id")
    parser.add_argument("--kingdom", choices=["Archaea", "Bacteria", "Mitochondria", "Viruses"],
                        help="Kingdom of the sample (default: Bacteria)", default="Bacteria"),
    parser.add_argument("--organism",
                        help="organism name"),
    parser.add_argument("--strain",
                        help="strain name"),
    parser.add_argument("--template",
                        help="NCBI template file"),
    parser.add_argument("--gcode", choices=GENETIC_CODES, type=int,
                        help="Genetic codes used, see "
                             "'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c'"
                             "for more information (default: 11)", default=11),
    parser = add_blast_args(parser)
    parser = add_workflow_args(parser)

    return parser


def add_all_args(parser):

    parser.add_argument("-p", "--prefix", metavar="FILE", required=True,
        help="Sample name.")
    parser.add_argument("--taxon", choices=["mitochondrion", "plasmid", "plastid", "viruses"], default="mitochondrion",
        help="Choose taxon for assembly species,default=plasmid.")
    parser.add_argument("-pl", "--platform", choices=["illumina", "mgi"], default="mgi",
        help="Select sequencing platform, default=mgi")
    parser.add_argument("-r1", "--reads1", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--reads2", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r2.")
    parser.add_argument("-s", "--sequencer", choices=["PromethION", "GridION", "RSII", "Sequel"], default="PromethION",
        help="Used sequencer, default=PromethION")
    parser.add_argument("-r", "--reads", metavar="FLIE", type=str, default="",
        help="Input the third-generation sequencing reads.")
    parser.add_argument("-d", "--database", metavar="FILE", type=str, default=None,
        help="Input the reference database path, default=None.")
    parser.add_argument("-pro", "--project", metavar="STR", type=str, default=None,
        help="Input project name, default=None.")
    parser.add_argument("-pid", "--projectid", metavar="STR", type=str, default=None,
        help="Input project id, default=None.")
    parser.add_argument("--trim", metavar="INT", type=int, default=5,
        help="Set trim length, default=5")
    parser.add_argument("--gc", metavar="INT", type=int, default=0,
        help="Set the maximum GC content retained, default=0")
    parser.add_argument("--base", metavar="STR", type=str, default="all",
        help="Set the number of reserved bases, default=all")
    parser.add_argument("--score", metavar="INT", type=int, default=0,
        help="Match score, default=0")
    parser.add_argument("--gcode", choices=GENETIC_CODES, type=int,
                        help="Genetic codes used, see "
                             "'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c'"
                             "for more information (default: 11)", default=11)
    parser.add_argument("-t", "--thread", type=int, default=4,
        help="Set the number of threads to run.")
    parser = add_workflow_args(parser)   

    return parser
