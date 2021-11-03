# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import os

report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"
configfile: os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"config/default_config.yaml")


rule cluster:
    input:
        "genecatalog/gene_catalog.faa",
        "genecatalog/clustering/orf2gene.h5"

rule subcluster:
    input:
        "genecatalog/clustering/orf2gene.h5",
        expand("genecatalog/subcluster/gc{id}.fasta",id=config['subclusterids']),
        expand("genecatalog/clustering/gene2gc{id}.h5",id=config['subclusterids'])



# add scripts
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"scripts"))
import utils






include: "rules/clustering.smk"



## add default resources
for r in workflow.rules:
    if not "mem" in r.resources:
        r.resources["mem"]=config["mem"]["default"]
    if not "time" in r.resources:
        r.resources["time"]=config["runtime"]["default"]
        
        
    # convert to new units
    r.resources["mem_mb"] = r.resources["mem"] * 1000
    r.resources["time_min"] = r.resources["time"] * 60

