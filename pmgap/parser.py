#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["add_ngs_qc_args", "add_all_args",
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
    parser.add_argument("--thread", type=int, default=4,
        help="Set the number of threads to run.")
    parser = add_workflow_args(parser)

    return parser


def add_all_args(parser):

    parser = add_ngs_qc_args(parser)

    return parser
