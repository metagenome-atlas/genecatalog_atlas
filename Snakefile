# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import os

report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"


rule cluster:
    input:
        "genecatalog/gene_catalog.faa",
        "genecatalog/clustering/orf2gene.tsv.gz",
        #"genecatalog/counts/median_coverage.tsv.gz",
rule annotate:
    input:
        "genecatalog/annotations/eggNog.tsv.gz"

# add scripts
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"scripts"))
import utils


configfile: os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"config/default_config.yaml")


DBDIR = os.path.realpath(config["database_dir"])
EGGNOG_DIR = os.path.join(DBDIR,'EggNOGV2')


include: "rules/clustering.smk"
include: "rules/eggNOG.smk"
#include: "quantify.smk"
