#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
[general]
project=一种xxxx细菌基因组测序
id=BJXWZ-201805002
sample=979
sequencer=GridION
author=张兴国
reviewer=XX
lib_num=1
cell_num=1
species=
strain=
"""
import sys
import json
import argparse
import logging
import os.path
import shutil
from datetime import datetime

from jinja2 import Template
from docxtpl import DocxTemplate

try:
    reload(sys)
    sys.setdefaultencoding( "utf-8" )
except:
    pass
try:
    from ConfigParser import ConfigParser
except:
    from configparser import ConfigParser


LOG = logging.getLogger(__name__)

__version__ = "1.2.3"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def check_path(path):

    path = os.path.abspath(path)

    if not os.path.exists(path):
        msg = "File not found '{path}'".format(**locals())
        LOG.error(msg)
        raise Exception(msg)

    return path


def check_paths(obj):
    """
    check the existence of paths
    :param obj:
    :return: abs paths
    """

    if isinstance(obj, list):
        r = []
        for path in obj:
            r.append(check_path(path))

        return r
    else:
        return check_path(obj)


def mkdir(d):
    """
    from FALCON_KIT
    :param d:
    :return:
    """
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        LOG.debug('mkdir {!r}'.format(d))
        os.makedirs(d)
    else:
        LOG.debug('mkdir {!r}, {!r} exist'.format(d, d))

    return d


def read_config(cfg):
    """
    read config fron ini
    :param cfg:
    :return:
    """
    check_paths(cfg)

    r = {}
    config = ConfigParser()
    config.read(cfg)

    for section in config.sections():
        r[section] = {}

        for option in config.options(section):
            value = config.get(section, option).strip()
            if type(value) == type(b''):
                value = value.decode("utf-8")
            r[section][option] = value

    return r


def read_tsv(file):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        r.append(line.split("\t"))

    return r


def read_table_reads(file):

    return read_tsv(file)


def read_table_contig(file):

    return read_tsv(file)[0]


def read_table_base(file):

    return read_tsv(file)


def read_table_gc_depth(gc_depth, fasta):

    r = []
    data = {}

    for line in open(fasta):
        line = line.strip()
        if line.startswith(">"):
            records = line.split()
            name = records[0][1:]
            for i in records:
                if i.startswith("[topology"):
                    data[name] = i[10:-1]
                else:
                    pass

    total_length = 0
    total_depth = 0
    for i in read_tsv(gc_depth):
        total_length += int(i[1])
        total_depth += int(i[1])*float(i[3])
        i[3] = "{:.4f}".format(float(i[3]))
        r.append(i+[data[i[0]]])
    depth = total_depth/total_length

    return r, depth


def read_table_codon(file):

    return read_tsv(file)


def read_table_structure(file):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        if "CDS" not in line:
            continue
        r.append(line.split("\t"))

    return r


def read_table_annotate(file):

    return read_tsv(file)


def run_report(cfg, fasta, jsons,
               table_reads, table_contig, table_base, table_length_gc, table_codon, table_structure, table_function,
               figure_base_content, figure_base_gc, figure_base_quality, figure_depth, figure_protein,
               tpl_html, out_dir
               ):

    out_dir = mkdir(out_dir)
    now = datetime.now()

    table_data = read_table_reads(table_reads)
    table_length_gc, depth = read_table_gc_depth(table_length_gc, fasta)

    r = {
        "author": "",
        "reviewer": "",
        "year": now.year,
        "month": now.month,
        "day": now.day,
        "project": "",
        "id": "",
        "name": "",
        "strain": "",
        "sequencer": "",
        "raw_bases": table_data[0][2],
        "cell_num": 1,
        "clean_bases": table_data[0][4],
        "cds_num": read_table_annotate(table_function)[-1][1],
        "depth": "{:,.4f}".format(depth),
        "table_data": table_data[0],
        "table_contig": read_table_contig(table_contig),
        "table_base": read_table_base(table_base),
        "table_gc_depth": table_length_gc,
        "table_codon": read_table_codon(table_codon),
        "table_structure": read_table_structure(table_structure),
        "table_function": read_table_annotate(table_function),
        "software": {},
        "database": {}

    }

    r.update(read_config(cfg)["general"])

    for i in jsons:
        with open(i) as fh:
            if isinstance(fh, bytes):
                fh = fh.decode("utf-8")
            js = json.load(fh)
            for k in js:
                r[k].update(js[k])

    figures = [figure_base_content, figure_base_gc, figure_base_quality, figure_depth, figure_protein]

    # html_report
    for i in ["images", "static"]:
        temp = os.path.join(out_dir, i)
        if os.path.exists(temp):
            shutil.rmtree(temp)
        shutil.copytree(os.path.join(tpl_html, i), temp)

    for i in figures:
        shutil.copy(i, os.path.join(out_dir, "images/"))


    line = open(os.path.join(tpl_html, "report.html")).read()
    if isinstance(line, bytes):
        line = line.decode("utf-8")
    tpl = Template(line)

    with open(os.path.join(out_dir, "report.html"), "w") as fh:
        line = tpl.render(r)

        if type(line) == type(b''):
            line = line.decode("utf-8")
        fh.write(line)

    return r


def report(args):

    run_report(
        cfg=args.cfg,
        fasta=args.fasta,
        jsons=args.jsons,
        table_reads=args.reads,
        table_contig=args.contig,
        table_base=args.base,
        table_length_gc=args.length_gc,
        table_codon=args.codon,
        table_structure=args.structure,
        table_function=args.function,
        figure_base_content=args.base_content,
        figure_base_gc=args.base_gc,
        figure_base_quality=args.base_quality,
        figure_depth=args.depth,
        figure_protein=args.protein,
        tpl_html=args.html,
        out_dir=args.out
    )


def add_report_args(parser):

    parser.add_argument("cfg", help="cfg of sample")
    parser.add_argument("--fasta", required=True, help="genomic.fasta")
    parser.add_argument("--out", required=True, help="")
    table_group = parser.add_argument_group(title="Tables", )
    table_group.add_argument("--reads", required=True, help="")
    table_group.add_argument("--contig", required=True, help="")
    table_group.add_argument("--base", required=True, help="")
    table_group.add_argument("--codon", required=True, help="")
    table_group.add_argument("--length_gc", required=True, help="")
    table_group.add_argument("--structure", required=True, help="")
    table_group.add_argument("--function", required=True, help="")

    json_group = parser.add_argument_group(title="Json", )
    json_group.add_argument("--jsons", nargs='+', required=True, help="")

    figure_group = parser.add_argument_group(title="Figure", )
    figure_group.add_argument("--base_content", required=True, help="")
    figure_group.add_argument("--base_gc", required=True, help="")
    figure_group.add_argument("--base_quality", required=True, help="")
    figure_group.add_argument("--depth", required=True, help="")
    figure_group.add_argument("--protein", required=True, help="")

    template_group = parser.add_argument_group(title="Template", )
    template_group.add_argument("--html", required=True, help="")

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


version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_report_args(parser)
    args = parser.parse_args()

    report(args)


if __name__ == "__main__":
    main()
