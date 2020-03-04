import os

localrules: input_genes
rule input_genes:
    input:
        faa= os.path.abspath(config['input_faa'])
    output:
        faa= temp("genecatalog/input.faa")
    shell:
        "ln -s {input} {output}"





rule createdb:
    input:
        "genecatalog/{catalogname}.faa"
    output:
        temp(directory("genecatalog/{catalogname}_mmseqdb"))
    threads:
        1
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/make_db/{catalogname}.log"
    benchmark:
        "logs/benchmarks/createdb/{catalogname}.tsv"
    shell:
        "mkdir {output} 2> {log} ; "
        "mmseqs createdb {input} {output}/db >> {log} 2>> {log} "




rule cluster_genes:
    input:
        db="genecatalog/input_mmseqdb"
    output:
        clusterdb = temp(directory("genecatalog/clustering/mmseqs"))
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/clustering/cluster_proteins.log"
    threads:
        config["threads"]
    params:
        tmpdir= os.path.join(config['tmpdir'],"mmseqs"),
        clustermethod = 'linclust' if config['clustermethod']=='linclust' else 'cluster',
        coverage=config['coverage'], #0.8,
        minid=config['minid'], # 0.00
        extra=config['extra'],
    shell:
        """
            mkdir -p {params.tmpdir} {output} 2>> {log}

            mmseqs {params.clustermethod} -c {params.coverage} \
            --min-seq-id {params.minid} {params.extra} \
            --threads {threads} {input.db}/db {output.clusterdb}/db {params.tmpdir}  >> {log} 2>> {log}

            rm -fr  {params.tmpdir} 2>> {log}
        """


rule get_rep_proteins:
    input:
        db= rules.cluster_genes.input.db,
        clusterdb = rules.cluster_genes.output.clusterdb,
    output:
        cluster_attribution = temp("genecatalog/orf2gene_oldnames.tsv"),
        rep_seqs_db = temp(directory("genecatalog/protein_catalog")),
        rep_seqs = temp("genecatalog/representatives_of_clusters.fasta")
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/clustering/get_rep_proteins.log"
    benchmark:
        "logs/benchmarks/get_rep_proteins.tsv"
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    shell:
        """
        mmseqs createtsv {input.db}/db {input.db}/db {input.clusterdb}/db {output.cluster_attribution}  > {log} 2>> {log}

        mkdir {output.rep_seqs_db} 2>> {log}

        mmseqs result2repseq {input.db}/db {input.clusterdb}/db {output.rep_seqs_db}/db  >> {log} 2>> {log}

        mmseqs result2flat {input.db}/db {input.db}/db {output.rep_seqs_db}/db {output.rep_seqs}  >> {log} 2>> {log}

        """



rule rename_gene_catalog:
    input:
        faa= "genecatalog/representatives_of_clusters.fasta",
        cluster_attribution = "genecatalog/orf2gene_oldnames.tsv",
    output:
        faa= "genecatalog/gene_catalog.faa",
        cluster_attribution = "genecatalog/clustering/orf2gene.tsv.gz",
    shadow: "minimal"
    benchmark:
        "logs/benchmarks/rename_catalog.tsv"
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    log:
        "logs/rename_catalog.log"
    params:
        prefix='Gene'
    script:
        "../scripts/rename_catalog.py"






def get_subcluster_id(wildcards):

    id= int(wildcards.id)

    assert (id>0) & (id<100), f"id should be an integer in [0,100], got {wildcards.id}"

    id = id/100

    if id >= float(config['minid']):
        logger.error("Id for gene subclustering should be lower than that the gene catalog"
                     f" {id} is not smaller than {config['minid']}"
                     )
        exit(1)
    return id

rule subcluster_genes:
    input:
        db=ancient("genecatalog/gene_catalog_mmseqdb"),
        faa="genecatalog/gene_catalog.faa" # used to update if genecatalog updates
    output:
        clusterdb = temp(directory("genecatalog/subcluster/mmseqs{id}")),
        tmpdir= temp(directory(os.path.join(config['tmpdir'],"subcluster{id}"))),
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/subcluster/cluster_{id}.log"
    threads:
        config["threads"]
    params:
        clustermethod = 'linclust' if config['clustermethod']=='linclust' else 'cluster',
        coverage=config['coverage'],
        minid= get_subcluster_id,
        extra=config['extra'],
    shell:
        """
            mkdir {output.tmpdir} 2> {log}
            mmseqs {params.clustermethod} -c {params.coverage} \
            --min-seq-id {params.minid} {params.extra} \
            --threads {threads} {input.db}/db {output.clusterdb}/db {output.tmpdir}  >>  {log} 2>> {log}
        """


rule get_rep_subclusters:
    input:
        db=ancient(rules.subcluster_genes.input.db),
        clusterdb = rules.subcluster_genes.output.clusterdb,
    output:
        cluster_attribution = temp("genecatalog/subcluster/gene2gc{id}_oldnames.tsv"),
        rep_seqs_db = temp(directory("genecatalog/subcluster/rep_gc{id}")),
        rep_seqs = temp("genecatalog/subcluster/representatives_gc{id}.fasta")
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/subcluster/get_rep_proteins_{id}.log"
    threads:
        config.get("threads", 1)
    shell:
        """
        mmseqs createtsv {input.db}/db {input.db}/db {input.clusterdb}/db {output.cluster_attribution}  > {log} 2>> {log}

        mkdir {output.rep_seqs_db} 2>> {log}

        mmseqs result2repseq {input.db}/db {input.clusterdb}/db {output.rep_seqs_db}/db  >> {log} 2>> {log}

        mmseqs result2flat {input.db}/db {input.db}/db {output.rep_seqs_db}/db {output.rep_seqs}  >> {log} 2>> {log}

        """
