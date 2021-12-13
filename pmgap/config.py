
import os.path
from collections import OrderedDict

ROOT = "/Work/pipeline/NMGAP/v1.3.0/"
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
}

# DATABASE CONFIGURE
NT_DATABASE = "/Work/database/nt/20210922/nt"


