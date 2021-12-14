#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from collections import OrderedDict
from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)
__version__ = "1.1.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_depth(file):

    data = OrderedDict()

    for line in read_tsv(file, '\t'):
        if line[0] not in data:
            data[line[0]] = []
        data[line[0]].append(int(line[2]))


    return data


def stat_coverage(depth_dict, prefix):

    fo = open("%s.window_depth.tsv" % prefix, "w")
    fo.write("#Contig\tSequencing depth (X)\n")

    for i in depth_dict:
        depth = depth_dict[i]
        fo.write("{0}\t{1:,.2f}\n".format(i, sum(depth)*1.0/len(depth)))

    fo.close()


def get_depth(file, prefix, window=100):

    data = read_depth(file)
    stat_coverage(data, prefix)
    r = {}

    for i in data:
        seqlen = len(data[i])
        if seqlen <= window:
            LOG.info("Sequence %s is too short" % i)
            continue

        ends = []
        depths = []
        for j in range(window, seqlen, window):
            ends.append(j)
            depths.append(sum(data[i][j-window:j+1])/window)
        if len(ends) <= 10:
            continue
        r[i] = [ends, depths]

    return r


def plot_depth(data, prefix):

    #height=len(data)//2+len(data)%2
    height=len(data)
    plt.style.use('ggplot')

    if len(data)>1:
        plt.figure(figsize=(12, height*4.5))
    else:
        plt.figure(figsize=(12, 4.5))
    matplotlib.rcParams['ytick.direction'] = 'out'
    matplotlib.rcParams['xtick.direction'] = 'out'

    i = 0
    for seqid in data:
        i=i+1
        if len(data)>1:
            ax = plt.subplot(height, 1, i)
        else:
            ax = plt.subplot(1,1,1)
#        matplotlib.rcParams['xtick.direction'] = 'out'
        ax.spines['top'].set_visible(False) #去掉上边框
        ax.spines['bottom'].set_visible(False) #去掉下边框
        ax.spines['left'].set_visible(False) #去掉左边框
        ax.spines['right'].set_visible(False) #去掉右边框
        plt.plot(data[seqid][0], data[seqid][1], label=seqid, linewidth ='1',color="#cb416b")
#        plt.legend(loc='upper right',frameon=False)
#        plt.legend("boxoff")
        plt.tight_layout(pad=2.5, w_pad=2.5, h_pad=1.0)
        font = {'weight': 'bold','size': 12,}

        ymax = max(data[seqid][1])
        ymin = 0
        xmax = max(data[seqid][0])

        plt.ylim(ymin, ymax*1.2)
        plt.xlabel("Position", font)
        plt.ylabel("Depth", font)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
    plt.savefig("%s.depth.png" % prefix, dpi=700)
    plt.savefig("%s.depth.pdf" % prefix)

    return 0


def add_help_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
                        help='output file of samtools depth.')
    parser.add_argument('-w', '--window', metavar='INT', type=int, default=100,
                        help='Set the override window size(bp), default=100.')
    parser.add_argument('-p', '--prefix', metavar='STR', type=str, default='out',
                        help='Output prefix, default=out.')

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
        plot_depth_stat.py -- Draw sequencing coverage map
attention:
        plot_depth_stat.py <inputfile> -w <window_size>
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    data = get_depth(args.input, args.prefix, args.window)

    plot_depth(data, args.prefix)


if __name__ == "__main__":
    main()
