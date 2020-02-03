
CONDAENV='../envs'
tempdir="/tmp"



query = config['compare_genecatalog']["query"]

rule search_all:
    input:
        expand("Genecatalog/comparisons/{query}-{target}/filtered_results.tsv",
               query=query,
               target= [k for k in config['compare_genecatalog']['catalogs'].keys() if k!= query]
               )

rule createdb:
    input:
        lambda wildcards: config['compare_genecatalog']['catalogs'][wildcards.catalogname]
    output:
        directory("{catalogname}_mmseqdb")
    threads:
        1
    conda:
        "%s/mmseqs.yaml" % CONDAENV
    shell:
        " mkdir {output}; "
        "mmseqs createdb {input} {output}/db "


rule search:
    input:
        target='{target}_mmseqdb',
        query='{query}_mmseqdb',
    output:
        db=directory("Genecatalog/comparisons/{query}-{target}/mmseqdb"),
    params:
        s=0.9,
        tmpdir= temp(directory(os.path.join(tempdir,"mmseqs_search")))
    conda:
        "%s/mmseqs.yaml" % CONDAENV
    threads:
        16
    log:
        "logs/Genecatalog/compare/search_{query}_{target}.txt"
    shell:
        "mkdir {output} {params.tmpdir} 2> {log};  "
        "mmseqs createindex {input.target}/db {params.tmpdir} >> {log} ;"
        "mmseqs search -s {params.s} "
        " --threads {threads} "
        "{input.query}/db {input.target}/db {output.db}/db {params.tmpdir} >> {log}"

#targetID  alnScore  seqIdentity  eVal  qStart  qEnd  qLen  tStart  tEnd  tLen
rule createtsv:
    input:
        target='{target}_mmseqdb',
        query='{query}_mmseqdb',
        resultdb= "Genecatalog/comparisons/{query}-{target}/mmseqdb",
    output:
        "Genecatalog/comparisons/{query}-{target}/all_results.tsv",
    log:
        "logs/Genecatalog/compare/createtsv_{query}_{target}.log"
    threads:
        1
    shell:
        "  mmseqs createtsv {input.query}/db {input.target}/db {input.resultdb}/db "
        "{output} "
        "--threads {threads} "
        " > >(tee   {log})"

localrules: filter
rule filter:
    input:
        rules.createtsv.output[0]
    output:
        "Genecatalog/comparisons/{query}-{target}/filtered_results.tsv"
    params:
        id_treshold=0.9
    run:
        with open(input[0]) as fin, open(output[0],'w') as fout:
            for line in fin:
                if float(line.split('\t')[3]) > params.id_treshold:
                    fout.write(line)
