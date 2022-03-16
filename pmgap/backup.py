#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import logging
import argparse
import os.path

from pmgap.config import *
from thirdparty.dagflow import DAG, Task, do_dag
from pmgap import __author__, __email__, __version__


LOG = logging.getLogger(__name__)


README = """
├── 01_Data	测序结果目录
│   ├── {name}.qc.xls	数据统计表
│   ├── {name}.base_quality.pdf	二代数据的碱基质量图
│   ├── {name}.base_quality.png	二代数据的碱基质量图
│   ├── {name}.base_content.png	二代数据的碱基分布图
│   ├── {name}.base_content.pdf	二代数据的碱基分布图
│   ├── {name}.base_gc.png	二代数据的gc分布图
│   ├── {name}.base_gc.pdf	二代数据的gc分布图
├── 02_Assembly	组装结果目录
│   ├── {name}.genomic.fasta	基因组序列fasta文件
│   ├── {name}.contig_stat.tsv	基因组序列长度统计表
│   ├── {name}.base_stat.tsv	基因组序列碱基分布统计表
│   ├── {name}.length_gc.xls	基因序列的长度和GC含量统计
│   ├── {name}.depth.png	基因组序列测序深度图
├── 03_Annotation	注释结果目录
│   ├── {name}.genomic.gff3	基因组注释gff3文件
│   ├── {name}.genomic.gb	基因组注释gb文件
│   ├── {name}.genomic.sqn	NCBI sequin文件
│   ├── {name}.RNA.fasta	基因组基因序列文件
│   ├── {name}.protein.fasta	基因组编码蛋白序列
│   ├── {name}.protein_length.png	编码蛋白序列长度分布图
│   ├── {name}.protein_length.pdf	编码蛋白序列长度分布图
│   ├── {name}.codon.tsv	密码子表
│   ├── {name}.structure_summary.tsv	基因组结构注释统计表
│   ├── {name}.function_summary.tsv	蛋白功能注释统计表
│   ├── {name}.merge.annotate.xls	蛋白功能注释汇总结果
│   ├── Function	功能注释结果目录
│         ├── {name}.prodigal.gff3	prodigal运行结果
│         ├── {name}.ipr.out	Interproscan运行结果
│         ├── {name}.Pfam.tsv	Pfam数据库注释结果
│         ├── {name}.TIGRFAMs.tsv	TIGRFAMs数据库注释
│         ├── {name}.GO.tsv	GO数据库注释结果
│         ├── {name}.WEGO.png	GO注释分类图
│         ├── {name}.WEGO.pdf	GO注释分类图
│         ├── {name}.other.tsv	其它数据库注释结果
│         ├── {name}.refseq.m6	refseq数据库blast比对结果
│         ├── {name}.refseq.tsv	refseq数据库注释结果
│         ├── {name}.species.tsv	refseq数据库的物种鉴定结果
│         ├── {name}.COG.xml	COG数据库blast比对结果
│         ├── {name}.COG.tsv	COG数据库注释结果
│         ├── {name}.COG.png	COG注释分类图
│         ├── {name}.COG.pdf	COG注释分类图
│         ├── {name}.KEGG.m6	KEGG数据库blast比对结果
│         ├── {name}.KEGG.tsv	KEGG数据库注释结果
│         ├── {name}.KEGG.png	KEGG注释分类图
│         ├── {name}.KEGG.pdf	KEGG注释分类图
│         ├── {name}.kegg_pathway.tsv	KEGG的pathway标记
│         ├── {name}.SwissProt.m6	SwissProt数据库blast比对结果
│         ├── {name}.SwissProt.tsv	SwissProt数据库注释结果
"""

def write_readme(result, name):

    with open('%s/readme.txt' % result, 'w') as fh:
        fh.write(README.format(name=name))


def create_md5sum_task(name, project, job_type, work_dir, out_dir):

    task = Task(
        id="md5sum",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 -q all.q",
        script="""
cd {outdir}
zip -r {project}_{name}.report.zip static/ images/ report.html
rm -rf static/ images/ report.html
cd {outdir}/result/
{bin}/makemd5 0* >result.md5
""".format(bin=BIN,
            name=name,
            project=project,
            outdir=out_dir
        )
    )

    return task


def backup(name, project, job_type, work_dir, out_dir):

    write_readme(os.path.join(out_dir, "result"), name)

    dag = DAG("backup")
    task = create_md5sum_task(name, project, job_type, work_dir, out_dir)

    dag.add_task(task)
    do_dag(dag, 1, 5)


if __name__ == "__main__":
    main()
