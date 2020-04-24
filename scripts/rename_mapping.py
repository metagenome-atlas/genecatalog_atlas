import sys
sys.stdout= open(snakemake.log[0],"w")
sys.stderr= open(snakemake.log[0],"a")

import pandas as pd

name_mapping= pd.read_csv(snakemake.input.name_mapping,index_col=0,sep='\t',squeeze=True)
assert type(name_mapping)==pd.Series

# read cluster mapping in chuncks
write_header=True
chuncknr= 0
for orf2gene in pd.read_csv(snakemake.input.cluster_mapping,
                            usecols=[0,1], #  clustermaping can have a tailing tab character leading to a
                           index_col=1, # the format is "{cluster}\t{orf}"
                           squeeze=True,
                           header=None,
                           sep='\t',
                           chunksize=1e7):

    assert type(orf2gene)==pd.Series
    orf2gene.name=snakemake.params.headers[1]
    orf2gene.index.name = snakemake.params.headers[0]

# map gene representative name to gene id, write to file with header only once

    orf2gene.map(name_mapping).to_csv(snakemake.output[0],sep='\t',header=write_header,mode='a')
    write_header=False

    chuncknr+=1
    print(f"processed chunck {chuncknr}")
