import os



localrules: filter_genes
rule filter_genes:
    input:
        faa=config['input_faa']
    output:
        faa= "genecatalog/all_genes/filtered_genes.faa",
    threads:
        1
    params:
        min_length=config['minlength']
    run:
        import pyfastx
        import io
        fa = pyfastx.Fasta(input[0])

        with open(output.faa,'w',buffering=io.DEFAULT_BUFFER_SIZE*1000) as out_faa:
            for gene in fa:
                if len(gene) >= params.min_length:
                    out_faa.write(gene.raw)

        os.remove(fa.file_name+'.fxi')



rule cluster_genes:
    input:
        faa= "genecatalog/all_genes/filtered_genes.faa"
    output:
        db=temp(directory("genecatalog/all_genes/predicted_genes")),
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


localrules: orf2gene
rule orf2gene:
    input:
        cluster_attribution = "genecatalog/orf2gene_oldnames.tsv",
    output:
        cluster_attribution = "genecatalog/clustering/orf2gene.tsv.gz",
    run:
        import pandas as pd
        # CLuterID    GeneID    empty third column
        orf2gene= pd.read_csv(input.cluster_attribution,index_col=1, header=None,sep='\t')

        protein_clusters_old_names= orf2gene[0].unique()

        map_names = dict(zip(protein_clusters_old_names,
                             utils.gen_names_for_range(len(protein_clusters_old_names),'Gene')))

        orf2gene['Gene'] = orf2gene[0].map(map_names)
        orf2gene.index.name='ORF'
        orf2gene['Gene'].to_csv(output.cluster_attribution,sep='\t',header=True)






localrules: rename_gene_catalog
rule rename_gene_catalog:
    input:
        orf2gene = "genecatalog/clustering/orf2gene.tsv.gz",
        representatives= "genecatalog/representatives_of_clusters.fasta"
    output:
        temp("genecatalog/representatives_of_clusters.fasta.fxi"),
        faa= "genecatalog/gene_catalog.faa"
    run:
        import pandas as pd
        import pyfastx
        import io
        from utils.fasta import str2multiline

        fa = pyfastx.Fasta(input.representatives)

        representatives= fa.keys()

        map_names= pd.read_csv(input.orf2gene,index_col=0,sep='\t').loc[representatives,'Gene']

        with open(output.faa,'w',buffering=io.DEFAULT_BUFFER_SIZE*1000) as outf:
            for gene in fa:
                new_name= map_names[gene.name]
                lines='\n'.join(str2multiline(gene.seq))
                outf.write(">{new_name} {gene.description}\n{lines}\n")
