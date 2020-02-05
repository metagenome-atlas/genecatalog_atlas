import os



rule cluster_genes:
    input:
        faa= config['input_faa']
    output:
        db=temp(directory("genecatalog/input_genes")),
        clusterdb = temp(directory("genecatalog/clustering/mmseqs"))
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/clustering/cluster_proteins.log"
    threads:
        config.get("threads", 1)
    params:
        tmpdir= os.path.join(config['tmpdir'],"mmseqs"),
        clustermethod = 'linclust' if config['clustermethod']=='linclust' else 'cluster',
        coverage=config['coverage'], #0.8,
        minid=config['minid'], # 0.00
        extra=config['extra'],
        clusterdb= lambda wc, output: os.path.join(output.clusterdb,'clusterdb'),
        db=lambda wc, output: os.path.join(output.db,'inputdb')
    shell:
        """
            mkdir -p {params.tmpdir} {output} 2>> {log}
            mmseqs createdb {input.faa} {params.db} &> {log}

            mmseqs {params.clustermethod} -c {params.coverage} \
            --min-seq-id {params.minid} {params.extra} \
            --threads {threads} {params.db} {params.clusterdb} {params.tmpdir}  >>  {log}

            rm -fr  {params.tmpdir} 2>> {log}
        """


rule get_rep_proteins:
    input:
        db= rules.cluster_genes.output.db,
        clusterdb = rules.cluster_genes.output.clusterdb,
    output:
        cluster_attribution = temp("genecatalog/orf2gene_oldnames.tsv"),
        rep_seqs_db = temp(directory("genecatalog/protein_catalog")),
        rep_seqs = temp("genecatalog/representatives_of_clusters.fasta")
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/clustering/get_rep_proteins.log"
    threads:
        config.get("threads", 1)
    params:
        clusterdb= lambda wc, input: os.path.join(input.clusterdb,'clusterdb'),
        db=lambda wc, input: os.path.join(input.db,'inputdb'),
        rep_seqs_db=lambda wc, output: os.path.join(output.rep_seqs_db,'db')
    shell:
        """
        mmseqs createtsv {params.db} {params.db} {params.clusterdb} {output.cluster_attribution}  > {log}

        mkdir {output.rep_seqs_db} 2>> {log}

        mmseqs result2repseq {params.db} {params.clusterdb} {params.rep_seqs_db}  >> {log}

        mmseqs result2flat {params.db} {params.db} {params.rep_seqs_db} {output.rep_seqs}  >> {log}

        """





localrules: rename_gene_catalog
rule rename_gene_catalog:
    input:
        faa= "genecatalog/representatives_of_clusters.fasta",
        cluster_attribution = "genecatalog/orf2gene_oldnames.tsv",
    output:
        faa= "genecatalog/gene_catalog.faa",
        cluster_attribution = "genecatalog/clustering/orf2gene.tsv.gz",
    shadow: "minimal"
    params:
        prefix='Gene'
    script:
        "../scripts/rename_catalog.py"




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
    shell:
        " mkdir {output}; "
        "mmseqs createdb {input} {output}/db "




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
            --threads {threads} {input.db}/db {output.clusterdb}/db {output.tmpdir}  &>  {log}
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
        mmseqs createtsv {input.db}/db {input.db}/db {input.clusterdb}/db {output.cluster_attribution}  > {log}

        mkdir {output.rep_seqs_db} 2>> {log}

        mmseqs result2repseq {input.db}/db {input.clusterdb}/db {output.rep_seqs_db}/db  >> {log}

        mmseqs result2flat {input.db}/db {input.db}/db {output.rep_seqs_db}/db {output.rep_seqs}  >> {log}

        """


rule subcluster:
    input:
        expand("genecatalog/subcluster/representatives_gc{id}.fasta",id=config['subclusterids'])
