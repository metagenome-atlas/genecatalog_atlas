import os

configfile: os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"config/default_config.yaml")




rule all:
    input:
        "genecatalog/compare/{comparison}_id{id}.tsv".format(comparison=
                                                      "_".join(config['compare_catalogs'].keys()),
                                                      id= config['compare_id']
                                                      )

rule createdb:
    input:
        list(config['compare_catalogs'].values())
    output:
        temp(directory("genecatalog/compare/{comparison}_mmseqdb"))
    threads:
        1
    conda:
        "envs/mmseqs.yaml"
    log:
        "logs/genecatalog/make_db/compare_{comparison}.log"
    benchmark:
        "logs/benchmarks/createdb/compare_{comparison}.tsv"
    shell:
        "mkdir {output} 2> {log} ; "
        "mmseqs createdb {input} {output}/db >> {log} 2>> {log} "



rule cluster_genes:
    input:
        db="genecatalog/compare/{comparison}_mmseqdb"
    output:
        clusterdb = temp(directory("genecatalog/compare/mmseqs_{comparison}"))
    conda:
        "envs/mmseqs.yaml"
    log:
        "logs/genecatalog/compare/cluster_proteins/{comparison}.log"
    threads:
        config["threads"]
    params:
        tmpdir= os.path.join(config['tmpdir'],"mmseqs"),
        clustermethod = 'linclust',
        coverage=config['coverage'], #0.8,
        minid=config['compare_id'], # 0.00
        extra=config['extra'],
    shell:
        """
            mkdir -p {params.tmpdir} {output} 2>> {log}

            mmseqs {params.clustermethod} -c {params.coverage} \
            --min-seq-id {params.minid} {params.extra} \
            --threads {threads} {input.db}/db {output.clusterdb}/db {params.tmpdir}  >> {log} 2>> {log}

            rm -fr  {params.tmpdir} 2>> {log}
        """



rule get_mapping_original:
    input:
        db= rules.cluster_genes.input.db,
        clusterdb = rules.cluster_genes.output.clusterdb,
    output:
        cluster_attribution = "genecatalog/compare/{comparison}.tsv",
    conda:
        "envs/mmseqs.yaml"
    log:
        "logs/genecatalog/compare/get_mapping_{comparison}.log"
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    shell:
        """
        mmseqs createtsv {input.db}/db {input.db}/db {input.clusterdb}/db {output.cluster_attribution}  > {log} 2>> {log}
        """






## add default resources
for r in workflow.rules:
    if not "mem" in r.resources:
        r.resources["mem"]=config["mem"]["default"]
    if not "time" in r.resources:
        r.resources["time"]=config["runtime"]["default"]
