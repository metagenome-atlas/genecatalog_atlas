JAVA_MEM_FRACTION=0.85

rule align_reads_to_Genecatalog:
    input:
        reads=get_quality_controlled_reads,
        fasta = "Genecatalog/gene_catalog.fna",
    output:
        sam = temp("Genecatalog/alignments/{sample}.sam")
    params:
        input = lambda wc, input : input_params_for_bbwrap( input.reads),
        maxsites = 4,
        ambiguous = 'all',
        minid = config['minid'],
        maxindel = 1 # default 16000 good for genome deletions but not necessarily for alignment to contigs
    shadow:
        "shallow"
    log:
        "logs/Genecatalog/alignment/{sample}_map.log"
    conda:
        "../envs/bbmap.yaml"
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    shell:
        """
        bbwrap.sh nodisk=t \
            local=t \
            ref={input.fasta} \
            {params.input} \
            trimreaddescriptions=t \
            out={output.sam} \
            threads={threads} \
            minid={params.minid} \
            mdtag=t \
            xstag=fs \
            nmtag=t \
            sam=1.3 \
            ambiguous={params.ambiguous} \
            secondary=t \
            saa=f \
            maxsites={params.maxsites} \
            -Xmx{resources.java_mem}G \
            2> {log}
        """


rule pileup_Genecatalog:
    input:
        sam = "Genecatalog/alignments/{sample}.sam",
        bam = "Genecatalog/alignments/{sample}.bam"
    output:
        covstats = temp("Genecatalog/alignments/{sample}_coverage.tsv"),
        basecov = temp("Genecatalog/alignments/{sample}_base_coverage.txt.gz"),
    params:
        pileup_secondary = 't' # a read maay map to different genes
    log:
        "logs/Genecatalog/alignment/{sample}_pileup.log"
    conda:
        "../envs/bbmap.yaml"
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    shell:
        """pileup.sh in={input.sam} \
               threads={threads} \
               -Xmx{resources.java_mem}G \
               covstats={output.covstats} \
               basecov={output.basecov} \
               secondary={params.pileup_secondary} \
                2> {log}
        """

localrules: combine_gene_coverages
rule combine_gene_coverages:
    input:
        covstats = expand("Genecatalog/alignments/{sample}_coverage.tsv",
            sample=SAMPLES)
    output:
        "Genecatalog/counts/median_coverage.tsv.gz",
        "Genecatalog/counts/Nmapped_reads.tsv.gz",
    run:

        import pandas as pd
        import os

        combined_cov={}
        combined_N_reads={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]
            data= pd.read_csv(cov_file,index_col=0,sep='\t')
            data.loc[data.Median_fold<0,'Median_fold']=0
            combined_cov[sample]= data.Median_fold
            combined_N_reads[sample] = data.Plus_reads+data.Minus_reads

        pd.DataFrame(combined_cov).to_csv(output[0],sep='\t',compression='gzip')
        pd.DataFrame(combined_N_reads).to_csv(output[1],sep='\t',compression='gzip')
