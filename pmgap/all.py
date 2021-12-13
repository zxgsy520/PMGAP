#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import os.path
import argparse
import logging
import shutil

from pmgap.config import *
from pmgap.common import read_config, mkdir, check_paths, rm
from pmgap.ngs_qc import run_ngs_qc
from pmgap.parser import add_all_args

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []

def run_all(reads1, reads2, prefix, job_type, work_dir, out_dir,
               concurrent=10, refresh=15, trim=3, thread=4):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    reads1 = check_paths(reads1)
    reads2 = check_paths(reads2) 

    work_dict = {
        "data": "01_Data",
    }
    
    result_dir = mkdir(os.path.join(out_dir, "result"))

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(result_dir, v))

    LOG.info("run qc")
    clean1, clean2, stat_qc, quality, content, gc, options = run_ngs_qc(
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

    with open(os.path.join(result_dir, "pmgap.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return 0


def _all(args):

    run_all(reads1=args.reads1,
            reads2=args.reads2,
            prefix=args.prefix,
            job_type=args.job_type,
            work_dir=args.work_dir,
            out_dir=args.out_dir,
            concurrent=args.concurrent,
            refresh=args.refresh, 
            trim=args.trim,
            thread=args.thread
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
