#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pmgap.parser import *
from pmgap.ngs_qc import ngs_qc

from pmgap.all import _all

from pmgap import __version__, __email__, __author__



def add_pmgap_parser(parser):

    subparsers = parser.add_subparsers(
        title='command',
        dest='commands')
    subparsers.required = True

    all_parser = subparsers.add_parser("all", help="all steps")
    all_parser = add_all_args(all_parser)
    all_parser.set_defaults(func=_all)

    ngs_qc_parser = subparsers.add_parser('ngs_qc', help="Quality control of second-generation data")
    ngs_qc_parser = add_ngs_qc_args(ngs_qc_parser)
    ngs_qc_parser.set_defaults(func=ngs_qc)

    return parser


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
URL:https://github.com/zxgsy520/PMGAP

describe:
Perform the assembly and annotation process of plasmids, mitochondria, plastids, and phages.

version: %s
contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    parser = add_pmgap_parser(parser)
    args = parser.parse_args()

    args.func(args)

    return parser.parse_args()



if __name__ == "__main__":
    main()
