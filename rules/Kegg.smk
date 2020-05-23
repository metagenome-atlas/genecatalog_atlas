

rule run_kofamscan:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        faa = "{folder}/{prefix}.faa",
    output:
        temp("{folder}/{prefix}.emapper.seed_orthologs"),
    params:
        data_dir = EGGNOG_DIR,
        prefix = "{folder}/{prefix}"
    resources:
        mem = config["mem"]["eggnog"],
        time = config["runtime"]["eggnog"]
    threads:
        config.get("threads_eggnog", config['threads'])
    benchmark:
        "{folder}/logs/benchmark/eggNOG_homology_search_diamond/{prefix}.log"
    conda:
        "../envs/eggNOG.yaml"
    log:
        "{folder}/logs/{prefix}/eggNOG_homology_search_diamond.log"
    shell:
        """
        emapper.py -m diamond --no_annot --no_file_comments \
            --data_dir {params.data_dir} --cpu {threads} -i {input.faa} \
            -o {params.prefix} --override 2> {log}
        """


def combine_kofam_annotations_input(wildcards):
    dir_for_subsets = checkpoints.gene_subsets.get(**wildcards).output[0]
    Subset_names,= glob_wildcards(os.path.join(dir_for_subsets, "{subset}.faa"))
    return expand("genecatalog/subsets/KOfamScan/{subset}.tsv",
                  subset=Subset_names)

rule combine_egg_nogg_annotations:
    input:
        combine_genecatalog_annotations_input
    output:
        "genecatalog/annotations/KO.tsv.gz"
    run:
        from utils.io import pandas_concat

        pandas_concat(input,output[0])
