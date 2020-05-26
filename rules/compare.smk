



rule createdb_compare:
    input:
        list(config['compare_catalogs'].values())
    output:
        temp(directory("genecatalog/compare_mmseqdb"))
    threads:
        1
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/make_db/compare.log"
    benchmark:
        "logs/benchmarks/createdb/compare.tsv"
    shell:
        "mkdir {output} 2> {log} ; "
        "mmseqs createdb {input} {output}/db >> {log} 2>> {log} "



rule compare_genes:
    input:
        db="genecatalog/compare_mmseqdb"
    output:
        clusterdb = temp(directory("genecatalog/clustering_compare/mmseqs"))
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/compare/cluster_proteins.log"
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



rule get_mapping_compare:
    input:
        db= rules.compare_genes.input.db,
        clusterdb = rules.compare_genes.output.clusterdb,
    output:
        cluster_attribution = "genecatalog/compare/"+ "_".join(config['compare_catalogs'].keys()) +".tsv",
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/compare/get_mapping.log"
    benchmark:
        "logs/benchmarks/get_mapping_compare.tsv"
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    shell:
        """
        mmseqs createtsv {input.db}/db {input.db}/db {input.clusterdb}/db {output.cluster_attribution}  > {log} 2>> {log}
        """


rule compare:
    input:
        rules.get_mapping_compare.output
