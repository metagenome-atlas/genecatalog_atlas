tempdir="/tmp"



query = config["query"]

rule search:
    input:
        expand("genecatalog/comparisons/{query}-{target}/filtered_results.tsv",
               query=query,
               target= [k for k in config['compare_catalogs'].keys() if k!= query]
               )


localrules: get_data_to_compare
rule get_data_to_compare:
    input:
        lambda wildcards: config['compare_catalogs'][wildcards.catalogname]
    output:
        temp("genecatalog/{catalogname}.faa")
    shell:
        "cp {input} {output}"

rule mmseq_search:
    input:
        target='genecatalog/{target}_mmseqdb',
        query='genecatalog/{query}_mmseqdb',
    output:
        db=directory("genecatalog/comparisons/{query}-{target}/mmseqdb"),
    params:
        s=0.9,
        tmpdir= temp(directory(os.path.join(tempdir,"mmseqs_search")))
    conda:
        "../envs/mmseqs.yaml"
    threads:
        config['threads']
    log:
        "logs/genecatalog/compare/search_{query}_{target}.txt"
    shell:
        "mkdir {output} {params.tmpdir} 2> {log};  "
        "mmseqs createindex {input.target}/db {params.tmpdir} >> {log} ;"
        "mmseqs search -s {params.s} "
        " --threads {threads} "
        "{input.query}/db {input.target}/db {output.db}/db {params.tmpdir} >> {log}"

#targetID  alnScore  seqIdentity  eVal  qStart  qEnd  qLen  tStart  tEnd  tLen
rule createtsv:
    input:
        target='genecatalog/{target}_mmseqdb',
        query='genecatalog/{query}_mmseqdb',
        resultdb= "genecatalog/comparisons/{query}-{target}/mmseqdb",
    output:
        "genecatalog/comparisons/{query}-{target}/all_results.tsv",
    log:
        "logs/genecatalog/compare/createtsv_{query}_{target}.log"
    threads:
        1
    shadow:
        "minimal"
    resources:
        mem=config['mem']['low'],
        time=config['runtime']['long']
    conda:
        "../envs/mmseqs.yaml"
    shell:
        "  mmseqs createtsv {input.query}/db {input.target}/db {input.resultdb}/db "
        "{output} "
        " > {log}"


rule filter:
    input:
        rules.createtsv.output[0]
    output:
        "genecatalog/comparisons/{query}-{target}/filtered_results.tsv"
    params:
        id_treshold=config['compare_id']
    threads:
        1
    resources:
        mem=config['mem']['low'],
    run:
        with open(input[0]) as fin, open(output[0],'w') as fout:
            for line in fin:
                if float(line.split('\t')[3]) > params.id_treshold:
                    fout.write(line)
