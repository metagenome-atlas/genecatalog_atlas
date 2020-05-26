import sys
sys.stdout= open(snakemake.log[0],"w")
sys.stderr= open(snakemake.log[0],"a")

import pandas as pd

name_mapping= pd.read_csv(snakemake.input.name_mapping,index_col=0,sep='\t',squeeze=True)
assert type(name_mapping)==pd.Series



with pd.HDFStore(snakemake.output[0],complevel=3, mode='w') as store:

    # read cluster mapping in chuncks
    chuncknr= 0
    for orf2gene in pd.read_csv(snakemake.input.cluster_mapping,
                                usecols=[0,1], #  clustermaping can have a tailing tab character leading to a
                               index_col=1, # the format is "{cluster}\t{orf}"
                               squeeze=True,
                               header=None,
                               sep='\t',
                               dtype={0:'category'},
                               chunksize=1e7
        ):


        orf2gene.cat.set_categories(name_mapping,
                                    inplace=True,
                                    rename=True,
                                    ordered=True)



        orf2gene.name=snakemake.params.headers[1]
        orf2gene.index.name = snakemake.params.headers[0]
        key= '/'.join(snakemake.params.headers)

    # map gene representative name to gene id, write to file with header only once
        store.append(key, orf2gene, format='table', data_columns=[orf2gene.name],
                     min_itemsize=50)

        chuncknr+=1
        print(f"processed chunck {chuncknr}")
