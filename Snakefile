# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


configfile: "config.yaml"
report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"


rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.



DBDIR = os.path.realpath(config["database_dir"])
EGGNOG_DIR = os.path.join(DBDIR,'EggNOGV2')
CONDAENV = "../envs"

include: "rules/other.smk"
