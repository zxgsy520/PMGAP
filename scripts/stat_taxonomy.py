#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep="\t"):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def stat_reads(file):
    
    n = 0
    for line in open(file):
        line = line.strip()
        
        if line.startswith(">"):
            n += 1

    return n


def stat_taxonomy(file, reads_number, name):

    total = reads_number
    map_number = 0
    species_dict = {}
    boundary_dict = {'Archaea':0, 'Bacteria':0 , 'Fungi':0, 'Viridiplantae':0, 'Alveolata':0, 'Metazoa':0 , 'Viruses':0, 'Viroids':0, 'Unclassified':0, 'Other':0}
    
    for line in read_tsv(file):
        category = "Other"
        classify = line[3].strip().split(';')

        map_number += 1
        if classify[0].strip() in boundary_dict:
             boundary_dict[classify[0].strip()] += 1
             category = classify[0].strip()

        elif classify[0].strip() == 'unclassified sequences':
             boundary_dict['Unclassified'] +=1
             category = classify[0].strip()

        elif classify[0].strip() == 'cellular organisms':
             if classify[1].strip() in boundary_dict:
                 boundary_dict[classify[1].strip()] += 1
                 category = classify[1].strip()

             elif classify[1].strip() == 'Eukaryota':
                 if classify[2].strip() in boundary_dict:
                     boundary_dict[classify[2].strip()] += 1
                     category = classify[2].strip()

                 elif classify[2].strip() == 'Opisthokonta':
                     if classify[3].strip() in boundary_dict:
                         boundary_dict[classify[3].strip()] += 1
                         category = classify[3].strip()

        else:
            boundary_dict['Other'] += 1

        species_name = '%s_%s' % (category, line[2])

        if species_name not in species_dict:
             species_dict[species_name] = 1
        else:
             species_dict[species_name] += 1
    
    out_bound = open('%s.species_classify.tsv' % name, 'w')
    out_bound.write('#Bacteria\tArchaea\tFungi\tViridiplantae\tAlveolata\tMetazoa\tViruses\tViroids\tUnclassified\tOther\n')
    out_bound.write('{:,}\t{:,}\t{:,}\t{:,}\t{:,}\t{:,}\t{:,}\t{:,}\t{:,}\t{:,}\n'.format(
        boundary_dict['Bacteria'],
        boundary_dict['Archaea'],
        boundary_dict['Fungi'],
        boundary_dict['Viridiplantae'],
        boundary_dict['Alveolata'],
        boundary_dict['Metazoa'],
        boundary_dict['Viruses'],
        boundary_dict['Viroids'],
        boundary_dict['Unclassified'],
        boundary_dict['Other'])
    )

    out_species = open('%s.stat_species.tsv' % name, 'w')
    out_species.write('Unknow\tUnmap\t{}\n'.format(total-map_number))

    for i in sorted(species_dict.items(),key = lambda x:x[1],reverse = True):
        species_name = i[0].split('_')
        out_species.write('{}\t{}\t{}\n'.format(species_name[0], species_name[1], i[1]))
    out_species.close()


def add_args(parser):

    parser.add_argument('-i', '--input', metavar='m6', type=str, required=True,
        help='Input the species annotation file.')
    parser.add_argument('-rn', '--reads_number', metavar='INT', type=int, default=10000,
        help='Number of reads participating in the comparison,default=10000.')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default= 'out',
        help='Output file prefix')
    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:

    stat_taxonomy.py -- Statistically annotated species information

attention:
    stat_taxonomy.py -i txt.species_annotation.txt -n name
    stat_taxonomy.py -i txt.species_annotation.txt -rn 1000000 -n out
''')
    args = add_args(parser).parse_args()

    stat_taxonomy(args.input, args.reads_number, args.name)


if __name__ == "__main__":
    main()
