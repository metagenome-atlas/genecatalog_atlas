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
        rep_seqs_db = temp(directory("genecatalog/protein_catalog")),
        rep_seqs = temp("genecatalog/representatives_of_clusters.fasta")
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/clustering/get_rep_proteins.log"
    benchmark:
        "logs/benchmarks/get_rep_proteins.tsv"
    resources:
        mem=config['mem']['low']
    threads:
        1
    shell:
        """

        mkdir {output.rep_seqs_db} 2>> {log}

        mmseqs result2repseq {input.db}/db {input.clusterdb}/db {output.rep_seqs_db}/db  >> {log} 2>> {log}

        mmseqs result2flat {input.db}/db {input.db}/db {output.rep_seqs_db}/db {output.rep_seqs}  >> {log} 2>> {log}

        """

rule get_mapping_original:
    input:
        db= rules.cluster_genes.input.db,
        clusterdb = rules.cluster_genes.output.clusterdb,
    output:
        cluster_attribution = temp("genecatalog/clustering/cluster_attribution.tsv"),
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/clustering/get_rep_proteins.log"
    benchmark:
        "logs/benchmarks/get_mapping_original.tsv"
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    shell:
        """
        mmseqs createtsv {input.db}/db {input.db}/db {input.clusterdb}/db {output.cluster_attribution}  > {log} 2>> {log}
        """



rule rename_gene_catalog:
    input:
        faa= "genecatalog/representatives_of_clusters.fasta",
        log= "logs/genecatalog/clustering/cluster_proteins.log"
    output:
        faa= "genecatalog/gene_catalog.faa",
        name_mapping = "genecatalog/clustering/renamed_genenames.tsv.gz",
    shadow: "minimal"
    benchmark:
        "logs/benchmarks/rename_catalog.tsv"
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    log:
        "logs/genecatalog/clustering/rename_catalog.log"
    params:
        prefix='Gene'
    script:
        "../scripts/rename_catalog.py"


rule rename_mapping:
    input:
        name_mapping = "genecatalog/clustering/renamed_genenames.tsv.gz",
        cluster_mapping ="genecatalog/clustering/cluster_attribution.tsv"
    output:
        "genecatalog/clustering/orf2gene.tsv.gz",
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    log:
        "logs/genecatalog/clustering/rename_mapping_clusters.log"
    shadow:
        "minimal"
    run:


        import pandas as pd

        name_mapping= pd.read_csv(input.name_mapping,index_col=0,sep='\t',squeeze=True)
        assert type(name_mapping)==pd.Series

        # read cluster mapping in chuncks
        write_header=True
        for orf2gene in pd.read_csv(input.cluster_mapping,
                                    usecols=[0,1], #  clustermaping can have a tailing tab character leading to a
                                   index_col=1, # the format is "{cluster}\t{orf}"
                                   squeeze=True,
                                   header=None,
                                   sep='\t',
                                   chunksize=1e7):

            assert type(orf2gene)==pd.Series
            orf2gene.name='Gene'
            orf2gene.index.name = 'ORF'

        # map gene representative name to gene id, write to file with header only once

            orf2gene.map(name_mapping).to_csv(output[0],sep='\t',header=write_header,mode='a')
            write_header=False


#### SUBCLUSTERING ####

def get_subcluster_id(wildcards):

    id= int(wildcards.id)

    assert (id>30) & (id<100), f"id should be an integer in [30,100], got {wildcards.id}"

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
        clusterdb = temp(directory("genecatalog/subcluster/GC{id}_mmseqs")),
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
        tmpdir= directory(os.path.join(config['tmpdir'],"GC{id}_subcluster"))
    shell:
        """
            mkdir -p {output} {params.tmpdir} 2> {log}
            mmseqs {params.clustermethod} -c {params.coverage} \
            --min-seq-id {params.minid} {params.extra} \
            --threads {threads} {input.db}/db {output.clusterdb}/db {params.tmpdir}  >>  {log} 2>> {log}
        """


rule get_rep_subclusters:
    input:
        db=ancient(rules.subcluster_genes.input.db),
        clusterdb = rules.subcluster_genes.output.clusterdb,
    output:
        rep_seqs_db = temp(directory("genecatalog/subcluster/GC{id}_rep")),
        rep_seqs = temp("genecatalog/subcluster/GC{id}_representatives.fasta")
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/subcluster/GC{id}_get_rep_proteins.log"
    threads:
        config.get("threads", 1)
    shell:
        """

        mkdir {output.rep_seqs_db} 2>> {log}

        mmseqs result2repseq {input.db}/db {input.clusterdb}/db {output.rep_seqs_db}/db  >> {log} 2>> {log}

        mmseqs result2flat {input.db}/db {input.db}/db {output.rep_seqs_db}/db {output.rep_seqs}  >> {log} 2>> {log}

        """

rule rename_subcluster_catalog:
    input:
        faa= "genecatalog/subcluster/GC{id}_representatives.fasta",
        log= "logs/genecatalog/subcluster/cluster_{id}.log"
    output:
        faa= "genecatalog/subcluster/gc{id}.fasta",
        name_mapping = temp("genecatalog/clustering/GC{id}_name_mapping.tsv"),
    shadow: "minimal"
    benchmark:
        "logs/benchmarks/GC{id}_rename_gene_clusters.tsv"
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    log:
        "logs/genecatalog/clustering/GC{id}_rename_gene_clusters.log"
    params:
        prefix='GC{id}_'
    script:
        "../scripts/rename_catalog.py"




rule get_subcluster_mapping_original:
    input:
        db= rules.subcluster_genes.input.db,
        clusterdb = rules.subcluster_genes.output.clusterdb,
    output:
        cluster_attribution = temp("genecatalog/clustering/GC{id}_cluster_attribution.tsv"),
    conda:
        "../envs/mmseqs.yaml"
    log:
        "logs/genecatalog/clustering/GC{id}_get_rep_proteins.log"
    benchmark:
        "logs/benchmarks/GC{id}_get_mapping_original.tsv"
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    shell:
        """
        mmseqs createtsv {input.db}/db {input.db}/db {input.clusterdb}/db {output.cluster_attribution}  > {log} 2>> {log}
        """

rule rename_subcluster_mapping:
    input:
        name_mapping = "genecatalog/clustering/GC{id}_name_mapping.tsv",
        cluster_mapping ="genecatalog/clustering/GC{id}_cluster_attribution.tsv"
    output:
        "genecatalog/clustering/gene2gc{id}.tsv.gz",
    resources:
        time=config['runtime']['long'],
        mem=config['mem']['low']
    threads:
        1
    log:
        "logs/genecatalog/clustering/GC{id}_rename_mapping_clusters.log"
    run:


        import pandas as pd

        name_mapping= pd.read_csv(input.name_mapping,index_col=0,sep='\t',squeeze=True)

        gene2gc= pd.read_csv(input.cluster_mapping,index_col=1,header=None,squeeze=True,sep='\t')

        assert type(name_mapping)==pd.Series
        assert type(gene2gc)==pd.Series

        gene2gc=pd.Series(index=gene2gc.index, data=name_mapping.loc[gene2gc.values].values,
                          name='GC{id}'.format(**wildcards)).sort_index()

        gene2gc.index.name='Gene'


        gene2gc.to_csv(output[0],sep='\t',header=True)
