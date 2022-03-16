PMGAP manual
===========
Table of Contents
-----------------

- [Quick usage](#quickusage)
- [Examples](#examples)


## <a name="quickusage"></a> Quick usage

```
./pmgap.py -h
usage: pmgap.py [-h] {all,ngs_qc,get_mito,assemble,annotate} ...

URL:https://github.com/zxgsy520/PMGAP

describe:
Perform the assembly and annotation process of plasmids, mitochondria, plastids, and phages.

version: 1.0.0
contact:  Xingguo Zhang <invicoun@foxmail.com>        

optional arguments:
  -h, --help            show this help message and exit

command:
  {all,ngs_qc,get_mito,assemble,annotate}
    all                 all steps              #运行整个流程
    ngs_qc              Quality control of second-generation data   #进行数据的质控和评估模块
    get_mito            Get mitochondrial chloroplast reads  #提取线粒体和叶绿体序列模块
    assemble            genome assembly #进行基因组的组装模块
    annotate            genome annotation #进行基因组的注释模块

```
```
./pmgap.py all -h
usage: pmgap.py all [-h] -p FILE
                    [--taxon {mitochondrion,plasmid,plastid,viruses}]
                    [-pl {illumina,mgi}] -r1 FILE [FILE ...] -r2 FILE
                    [FILE ...] [-s {PromethION,GridION,RSII,Sequel}] [-r FLIE]
                    [-d FILE] [-pro STR] [-pid STR] [--trim INT] [--gc INT]
                    [--base STR] [--score INT]
                    [--gcode {1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31}]
                    [-t THREAD] [--concurrent INT] [--refresh INT]
                    [--job_type {sge,local}] [--work_dir DIR] [--out_dir DIR]

optional arguments:
  -h, --help            show this help message and exit
  -p FILE, --prefix FILE  #输入样本的名称
                        Sample name.
  --taxon {mitochondrion,plasmid,plastid,viruses} #分析的物种类型
                        Choose taxon for assembly species,default=plasmid.
  -pl {illumina,mgi}, --platform {illumina,mgi} #分析数据使用的测序平台
                        Select sequencing platform, default=mgi
  -r1 FILE [FILE ...], --reads1 FILE [FILE ...] #输入测序数据的R1，允许通配符和多个文件
                        Second generation sequencing of reads r1.
  -r2 FILE [FILE ...], --reads2 FILE [FILE ...] #输入测序数据的R2，允许通配符和多个文件
                        Second generation sequencing of reads r2.
  -s {PromethION,GridION,RSII,Sequel}, --sequencer {PromethION,GridION,RSII,Sequel} #三代测序数据的平台
                        Used sequencer, default=PromethION
  -r FLIE, --reads FLIE #三代测序数据，非必须
                        Input the third-generation sequencing reads.
  -d FILE, --database FILE  #参考数据库，用于提取目标序列
                        Input the reference database path, default=None.
  -pro STR, --project STR #项目名称
                        Input project name, default=None.
  -pid STR, --projectid STR #项目编号
                        Input project id, default=None.
  --trim INT            Set trim length, default=5  #二代数据指控使用的修剪的碱基数目
  --gc INT              Set the maximum GC content retained, default=0  #对于没有比对上数据库的reads，保留的最高GC含量
  --base STR            Set the number of reserved bases, default=all #提取的碱基数目，数据量太大，可以用这个值保留固定的数据量进行后续组装
  --score INT           Match score, default=0 #比对的得分
  --gcode {1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31}  #基因组对应的密码子表
                        Genetic codes used, see 'https://www.ncbi.nlm.nih.gov/
                        Taxonomy/Utils/wprintgc.cgi?mode=c'for more
                        information (default: 11)
  -t THREAD, --thread THREAD  #运行的线程数目
                        Set the number of threads to run.

Workflow arguments:
  --concurrent INT      Maximum number of jobs concurrent (default: 10) #平行的任务数目
  --refresh INT         Refresh time of log in seconds (default: 30)  #检测任务状态的刷新时间
  --job_type {sge,local}
                        Jobs run on [sge, local] (default: local) #sge为使用sge任务调度系统，local为本地运行
  --work_dir DIR        Work directory (default: current directory) #流程的工作路径
  --out_dir DIR         Output directory (default: current directory) #流程的结果路径
```
## <a name="examples"></a> Examples
### single sample analysis
```
/Work/pipeline/PMGAP/v1.0.0/pmgap.py all \
  --reads1 raw_1.fq.gz \
  --reads2 raw_2.fq.gz \
  --prefix NC --taxon viruses --trim 5 --thread 10 --job_type sge \
  --platform mgi --database viral.genomic.fasta \
  --project 病毒基因组测序分析 --projectid NC2022010501 \
  --work_dir work  --out_dir out --gcode 1
```
### Multiple Sample Analysis
```
/Work/pipeline/PMGAP/v1.0.0/pmgap/pmgap_multi.py ngs.list --taxon viruses \   #ngs.list为二代数据路径文件
  --job_type sge --work_dir work --out_dir out --platform mgi \
  --gcode 1 --database viral.genomic.fasta \   #viral.genomic.fasta为参考基因组文件，可以不输入参考，将使用全部数据进行组装
  --project 病毒基因组测序分析 --projectid NC2022010501

ngs.list format：
NC1	/Project/data/NC1_raw_1.fq.gz	/Project/data/NC1_raw_2.fq.gz   #第一列为样本名称，第二列为二代测序的R1文件绝对路径，第三列为二代测序的R1文件绝对路径，
NC2	/Project/data/NC2_raw_1.fq.gz	/Project/data/NC2_raw_2.fq.gz
```
