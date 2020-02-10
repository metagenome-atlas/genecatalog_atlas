import pandas as pd
import utils
import sys
sys.stdout= open(snakemake.log[0],"w")
sys.stderr= open(snakemake.log[0],"a")




def rename_fasta(fasta_in,fasta_out,map_names):

    assert type(map_names)==dict

    with open(fasta_in) as fin, open(fasta_out,'w') as fout:
        for line in fin:
            if line[0]=='>':
                name=line[1:].split(maxsplit=1)[0]
                new_name=map_names.pop(name)
                line=f">{new_name}"

            fout.write(line)




if __name__ == '__main__':


    # Rename mapping
    #  ORF to gene

    # Noheaders CLuterID    GeneID    empty third column
    orf2gene= pd.read_csv(snakemake.input.cluster_attribution,index_col=1, header=None,sep='\t').iloc[:,0]


    representatives= orf2gene.unique()
    gene_names=utils.gen_names_for_range(len(representatives),snakemake.params.prefix)
    map_names = dict(zip(representatives,gene_names))

    orf2gene = orf2gene.map(map_names)
    orf2gene.index.name='ORF'
    orf2gene.name = 'Gene'


    duplicated_genes=[]
    if not orf2gene.index.is_unique:
        raise Exception("ORF names are not unique, remove duplicates see:\n",
              orf2gene.loc[orf2gene.index.duplicated(keep=False)])


        #orf2gene=orf2gene.loc[~orf2gene.index.duplicated()]


    orf2gene.to_csv(snakemake.output.cluster_attribution,sep='\t',header=True)

    del orf2gene

    # Rename representative sequence

    rename_fasta(snakemake.input.faa,snakemake.output.faa,map_names)
