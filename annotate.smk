import os



# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"
configfile: os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"config/default_config.yaml")



rule annotate:
    input:
        "genecatalog/annotations/eggNog.tsv.gz"

# add scripts
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"scripts"))
import utils




DBDIR = os.path.realpath(config["database_dir"])
EGGNOG_DIR = os.path.join(DBDIR,'EggNOGV2')



include: "rules/eggNOG.smk"



## add default resources
for r in workflow.rules:
    if not "mem" in r.resources:
        r.resources["mem"]=config["mem"]["default"]
    if not "time" in r.resources:
        r.resources["time"]=config["runtime"]["default"]
