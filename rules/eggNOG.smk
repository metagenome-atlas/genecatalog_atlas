
# Download
import hashlib
import os


# note: saving OG_fasta.tar.gz in order to not create secondary "success" file
FILES = {"eggnog.db": "7923d3bb7eca8e0e8f122be4b5ca6997",
         "eggnog_proteins.dmnd": "64fefa838833a6f3e220a06fb9d403cd"
         }

EGGNOG_HEADER = [
    "Query",
    "Target",
    "evalue",
    "score",
    "Taxonomy",
    "Protein_name",
    "GO_terms",
    "EC",
    "KO",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "tax_scope",
    "EggNog",
    "depricated_bestOG",
    "FunctionalCategory",
    "Description"
                ]


def get_eggnog_db_file():
    return ancient(expand("{path}/{files}",
                  path=EGGNOG_DIR,
                  files=["eggnog.db","eggnog_proteins.dmnd","checksum_checked"]
                  ))



def md5(fname):
    # https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    if not os.path.exists(fname):
        return None
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


localrules: download_eggNOG_files, verify_eggNOG_files


rule download_eggNOG_files:
    output:
        f"{EGGNOG_DIR}/eggnog.db",
        f"{EGGNOG_DIR}/eggnog_proteins.dmnd"
    threads:
        1
    conda:
        "../envs/eggNOG.yaml"
    shell:
        f"download_eggnog_data.py -yf --data_dir {EGGNOG_DIR} "

rule verify_eggNOG_files:
    input:
        rules.download_eggNOG_files.output
    output:
        touch(f"{EGGNOG_DIR}/checksum_checked")
    run:
        # validate the download
        for file in input:
            if not FILES[os.path.basename(file)] == md5(file):
                raise OSError(2, "Invalid checksum", file)
        # check if old eggNOG dir exists
        old_eggnogdir=EGGNOG_DIR.replace('V2','')
        if os.path.exists(old_eggnogdir):
            logger.info("The eggnog database form the olf v1 was found on your system."
                        f"You can savely remove this folder {old_eggnogdir}")

# RUN egg NOG
# # this rule specifies the more general eggNOG rules

# output with wildcards "{folder}/{prefix}.emapper.tsv"
folder= 'genecatalog/subsets/genes'

rule eggNOG_homology_search:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        faa = f"{folder}/{{prefix}}.faa",
    output:
        temp(f"{folder}/{{prefix}}.emapper.seed_orthologs"),
    params:
        data_dir = EGGNOG_DIR,
        prefix = f"{folder}/{{prefix}}"
    resources:
        mem = config["mem"]["eggnog"],
        time = config["runtime"]["eggnog"]
    threads:
        config.get("threads_eggnog", config['threads'])
    benchmark:
        "logs/benchmark/eggNOG/eggNOG_homology_search_diamond/{prefix}.tsv"
    conda:
        "../envs/eggNOG.yaml"
    log:
        "logs/eggNOG/eggNOG_homology_search_diamond/{prefix}.log"
    shell:
        """
        emapper.py -m diamond --no_annot --no_file_comments \
            --data_dir {params.data_dir} --cpu {threads} -i {input.faa} \
            -o {params.prefix} --override 2> {log}
        """



rule eggNOG_annotation:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        seed = rules.eggNOG_homology_search.output
    output:
        temp(f"{folder}/{{prefix}}.emapper.annotations")
    params:
        data_dir = EGGNOG_DIR,
        prefix = f"{folder}/{{prefix}}"
    threads:
        config.get("threads_eggnog", config['threads'])
    resources:
        mem=config['mem']['eggnog'],
        time = config["runtime"]["eggnog"]
    conda:
        "../envs/eggNOG.yaml"
    log:
        "logs/eggNOG/annotate_hits_table/{prefix}.log"
    benchmark:
        "logs/benchmark/eggNOG/annotate_hits_table/{prefix}.tsv"
    shell:
        """
        cp {params.data_dir}/eggnog.db /dev/shm/eggnog.db 2> {log}

        emapper.py --annotate_hits_table {input.seed} --no_file_comments \
            --override -o {params.prefix} --cpu {threads} --data_dir /dev/shm/ 2>> {log}
        """

localrules: add_eggNOG_header
rule add_eggNOG_header:
    input:
        f"{folder}/{{prefix}}.emapper.annotations"
    output:
        f"{folder}/{{prefix}}.eggNOG.tsv.gz"
    run:
        import pandas as pd

        D = pd.read_csv(input[0], header=None,sep='\t')
        D.columns = EGGNOG_HEADER
        D.to_csv(output[0],sep='\t')



localrules: gene_subsets
checkpoint gene_subsets:
    input:
        "genecatalog/gene_catalog.faa"
    output:
        directory(folder)
    params:
        subset_size=config['SubsetSize'],
    run:
        from utils import fasta
        fasta.split(input[0],params.subset_size,output[0])


def combine_genecatalog_annotations_input(wildcards):
    dir_for_subsets = checkpoints.gene_subsets.get(**wildcards).output[0]
    Subset_names,= glob_wildcards(os.path.join(dir_for_subsets, "{subset}.faa"))
    return expand("{folder}/{subset}.eggNOG.tsv.gz",
                  subset=Subset_names,folder=folder)

rule combine_egg_nogg_annotations:
    input:
        combine_genecatalog_annotations_input
    output:
        "genecatalog/annotations/eggNog.tsv.gz"
    run:
        from utils.io import pandas_concat

        pandas_concat(input,output[0])
