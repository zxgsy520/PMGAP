
import os.path
from collections import OrderedDict

ROOT = "/Work/pipeline/PMGAP/v1.0.0/"
SCRIPTS = os.path.join(ROOT, "scripts")
BIN = os.path.join(ROOT, "bin")
TEMPLATES = os.path.join(ROOT, "template")
TAXONOMY = os.path.join(TEMPLATES, "species.taxonomy")
QUEUE = "-q all.q,s01"

PYTHON_BIN = "/Work/pipeline/software/Base/miniconda3/bin/"

#ngs_qc
FASTP_BIN = "/Work/pipeline/software/Base/fastp/lastest/"
FASTQC_BIN = "/Work/pipeline/software/Base/fastqc/lastest/"
BLAST_BIN = "/Work/pipeline/software/Base/blast+/bin/"
#assembly
UNICYCLER_BIN = "/Work/pipeline/software/meta/unicycler/v0.4.8p/bin/"
BWA_BIN = "/Work/pipeline/software/Base/bwa/lastest/"
MUMMER_BIN = "/Work/pipeline/software/Base/mummer/lastest/bin/"
PRODIGAL_BIN = "/Work/pipeline/software/meta/prokka/lastest/bin/"
CIRCLATOR_BIN = "/Work/pipeline/software/meta/circlator/lastest/bin/"
SAMTOOLS_BIN = "/Work/pipeline/software/Base/samtools/samtools/"
MINIMAP2_BIN = "/Work/pipeline/software/Base/minimap2/"


SOFTWARE_VERSION = {
    "fastp":{
        "GETVER": "%s/fastp --version 2>&1 |grep 'fastp'" % FASTP_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.20.0",
    },
    "fastqc":{
        "GETVER": "%s/fastqc --version 2>&1 |grep 'FastQC'" % FASTQC_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.11.8",
    },
    "blastn":{
        "GETVER": "%s/blastn -version 2>&1| grep 'blastn'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+"
    },
    "unicycler": {
        "GETVER": "%s/unicycler --version 2>&1| grep 'Unicycler'" % UNICYCLER_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.4.9"
    },
    "bwa": {
        "GETVER": "%s/bwa 2>&1|grep -i '^Version:'" % BWA_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.7.17"
    },
    "samtools": {
        "GETVER": "%s/samtools 2>&1|grep -i '^Version:'" % SAMTOOLS_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.9"
    },
    "circlator": {
        "GETVER": "%s/circlator.py version" % CIRCLATOR_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.5.5"
    },
    "minimap2": {
        "GETVER": "%s/minimap2 --version" % MINIMAP2_BIN,
        "REGEXP": "\d+\.\d+\-r\d+",
        "MINVER": "2.11-r797"
    },
}

# DATABASE CONFIGURE
NT_DATABASE = "/Work/database/nt/20210922/nt"


