import pandas as pd
import utils


def rename_fasta(fasta_in,fasta_out,map_names):

    assert type(map_names)==dict

    with open(fasta_in) as fin, open(fasta_out,'w') as fout:
        for line in fin:
            if line[0]=='>':
                name=line[1:].split(maxsplit=1)[0]
                new_name=map_names.pop(name)
                line=line.replace(name,new_name)

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
    orf2gene.to_csv(snakemake.output.cluster_attribution,sep='\t',header=True)


    fasta_in=snakemake.input.faa
    fasta_out=snakemake.output.faa



    rename_fasta(fasta_in,fasta_out,map_names)
