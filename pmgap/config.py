
import os.path
from collections import OrderedDict

ROOT = "/Work/pipeline/PMGAP/v1.0.0/"
SCRIPTS = os.path.join(ROOT, "scripts")
BIN = os.path.join(ROOT, "bin")
TEMPLATES = os.path.join(ROOT, "template")
DATABASE = os.path.join(ROOT, "database")
TAXONOMY = os.path.join(DATABASE, "species.taxonomy.gz")
QUEUE = "-q all.q,s01"

GENETIC_CODES = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, ]
PYTHON_BIN = "/Work/pipeline/software/Base/miniconda3/bin/"

#ngs_qc
FASTP_BIN = "/Work/pipeline/software/Base/fastp/lastest/"
FASTQC_BIN = "/Work/pipeline/software/Base/fastqc/lastest/"
BLAST_BIN = "/Work/pipeline/software/Base/blast+/bin/"

SEQUENCER = OrderedDict([
    ("RSII", {
        "canu": "-pacobio-raw",
        "flye": "--plasmids --pacbio-raw",
        "minimap2": "-ax map-pb"}
     ),
    ("Sequel", {
        "canu": "-pacobio-raw",
        "flye": "--plasmids --pacbio-raw",
        "minimap2": "-ax map-pb"}
     ),
    ("GridION", {
        "canu": "-nanopore-raw",
        "flye": "--plasmids --nano-raw",
        "minimap2": "-ax map-ont"}
     ),
    ("PromethION", {
        "canu": "-nanopore-raw",
        "flye": "--plasmids --nano-raw",
        "minimap2": "-ax map-ont"}
     ),
    ])

#assembly
UNICYCLER_BIN = "/Work/pipeline/software/meta/unicycler/v0.4.8p/bin/"
BWA_BIN = "/Work/pipeline/software/Base/bwa/lastest/"
MUMMER_BIN = "/Work/pipeline/software/Base/mummer/lastest/bin/"
PRODIGAL_BIN = "/Work/pipeline/software/meta/prokka/lastest/bin/"
CIRCLATOR_BIN = "/Work/pipeline/software/meta/circlator/lastest/bin/"
SAMTOOLS_BIN = "/Work/pipeline/software/Base/samtools/samtools/"
MINIMAP2_BIN = "/Work/pipeline/software/Base/minimap2/"

#ANNOTATION
#PRODIGAL_BIN = "/Work/pipeline/software/meta/prokka/lastest/bin/"
GENEMARKS_BIN = "/Work/pipeline/software/meta/GeneMarkS/v4.25/"
VGAS_BIN = "/Work/pipeline/software/meta/vgas/v2020.07/bin/"
TBL2ASN_BIN = "/Work/pipeline/software/meta/tbl2asn/"
RSCRIPT = "Rscript"
DIAMOND_BIN = "/Work/pipeline/software/Base/diamond/v2.0.3/"
INTERPROSCAN_BIN = "/Work/pipeline/software/RNAseq/interproscan-5.51-85.0/"
ABRICATE_BIN = "/Work/pipeline/software/meta/abricate/v1.0.1/bin/" 

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
    "prodigal": {
        "GETVER": "%s/prodigal -v 2>&1| grep 'Prodigal'" % PRODIGAL_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.6.3"
    },
    "blastp": {
        "GETVER": "%s/blastp -version 2>&1| grep 'blastp'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+"
    },
    "rpsblast": {
        "GETVER": "%s/rpsblast -version 2>&1| grep 'rpsblast'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+"
    },
    "interproscan": {
        "GETVER": "%s/interproscan.sh -version 2>&1| grep 'version'" % INTERPROSCAN_BIN,
        "REGEXP": "\d+\.\d+\-\d+\.\d+",
        "MINVER": "5.30-69.0"
    },
    "tbl2asn": {
        "GETVER": "%s/tbl2asn --help | grep 'tbl2asn'" % TBL2ASN_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "25.6"
    },
}

# DATABASE CONFIGURE
NT_DATABASE = "/Work/database/nt/20210922/nt"
REFSEQ = "/Work/database/Refseq/202109/"
REFSEQ_TAXON = {
    "plasmid": os.path.join(REFSEQ, "protein/split/bacteria.dmnd"),
    "phage": os.path.join(REFSEQ, "protein/split/bacteria.dmnd"),
    "viruses": os.path.join(REFSEQ, "protein/split/viruses.dmnd"),
    "mitochondrion":  os.path.join(REFSEQ, "genomes/mitochondrion/mitochondrion.genomic.fasta")
}

KEGG = "/Work/database/kegg/2021/"
KEGG_KEG = os.path.join(KEGG, "ko00001.keg")
KEGG_PATH = os.path.join(KEGG, "ko2pathway.tsv")
KEGG_KO = os.path.join(KEGG, "microbe/microbe.kegg.pep2ko.txt")
KEGG_TAXON = {
    "plasmid": os.path.join(KEGG, "microbe/microbe.dmnd"),
    "phage": os.path.join(KEGG, "microbe/microbe.dmnd"),
    "viruses": os.path.join(KEGG, "microbe/microbe.dmnd"),
    "mitochondrion":""
}

COG = "/Work/database/COG/2020/"
COG_TAXON = {
    "plasmid": os.path.join(COG, "cog.dmnd"),
    "phage": os.path.join(COG, "cog.dmnd"),
    "viruses": os.path.join(COG, "cog.dmnd"),
    "mitochondrion": ""

}
COG_NAME = os.path.join(COG, "cog-20.def.tab")
COG_FUNC = os.path.join(COG, "fun-20.tab")

GO = "/Work/database/GO/20180828/"
TIGRFAMS = "/Work/database/TIGRFAMs/"
PFAM = "/Work/database/Pfam/v34.0/"
SWISSPROT_DB = "/Work/database/SwissProt/202111/swissprot.dmnd"
