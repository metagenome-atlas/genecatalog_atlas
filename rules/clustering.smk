import os



# localrules: concat_genes
# rule concat_genes:
#     input:
#         faa= expand("{sample}/annotation/predicted_genes/{sample}.faa", sample=SAMPLES),
#         fna= expand("{sample}/annotation/predicted_genes/{sample}.fna", sample=SAMPLES)
#     output:
#         faa=  temp("genecatalog/all_genes_unfiltered.faa"),
#         fna = temp("genecatalog/all_genes_unfiltered.fna"),
#     run:
#         from utils.io import cat_files
#         cat_files(input.faa,output.faa)
#         cat_files(input.fna,output.fna)
#
#



localrules: filter_genes
rule filter_genes:
    input:
        fna=config['input_fna'],
        faa=config['input_faa']
    output:
        fna= "genecatalog/all_genes/predicted_genes.fna",
        faa= "genecatalog/all_genes/predicted_genes.faa",
    threads:
        1
    params:
        min_length=config['minlength_nt']
    run:
        from Bio import SeqIO
        faa = SeqIO.parse(input.faa,'fasta')
        fna = SeqIO.parse(input.fna,'fasta')

        with open(output.faa,'w') as out_faa, open(output.fna,'w') as out_fna:

            for gene in fna:
                protein = next(faa)

                if len(gene) >= params.min_length:
                    SeqIO.write(gene,out_fna,'fasta')
                    SeqIO.write(protein,out_faa,'fasta')



rule cluster_genes:
    input:
        faa= "genecatalog/all_genes/predicted_genes.faa"
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
            --threads {threads} {params.db} {params.clusterdb} {params.tmpdir}  &>>  {log}

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
        mmseqs createtsv {params.db} {params.db} {params.clusterdb} {output.cluster_attribution}  &> {log}

        mkdir {output.rep_seqs_db} 2>> {log}

        mmseqs result2repseq {params.db} {params.clusterdb} {params.rep_seqs_db}  &>> {log}

        mmseqs result2flat {params.db} {params.db} {params.rep_seqs_db} {output.rep_seqs}  &>> {log}

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
        fna = "genecatalog/all_genes/predicted_genes.fna",
        faa= "genecatalog/all_genes/predicted_genes.faa",
        orf2gene = "genecatalog/clustering/orf2gene.tsv.gz",
        representatives= "genecatalog/representatives_of_clusters.fasta"
    output:
        fna= "genecatalog/gene_catalog.fna",
        faa= "genecatalog/gene_catalog.faa",
    run:
        import pandas as pd
        from Bio import SeqIO

        representatives= []
        with open(input.representatives) as fasta:
            for line in fasta:
                if line[0]=='>': representatives.append(line[1:].split()[0])

        map_names= pd.read_csv(input.orf2gene,index_col=0,sep='\t').loc[representatives,'Gene']

        # rename fna
        faa_parser = SeqIO.parse(input.faa,'fasta')
        fna_parser = SeqIO.parse(input.fna,'fasta')

        with open(output.fna,'w') as fna, open(output.faa,'w') as faa :
            for gene in fna_parser:
                protein = next(faa_parser)
                if gene.name in map_names.index:
                    gene.id = map_names[gene.name]
                    protein.id = map_names[protein.name]

                    SeqIO.write(gene,fna,'fasta')
                    SeqIO.write(protein,faa,'fasta')
