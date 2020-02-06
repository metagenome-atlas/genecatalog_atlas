import pandas as pd
import utils
import sys
sys.stdout= open(snakemake.log[0],"w")
sys.stderr= open(snakemake.log[0],"a")


from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of name as string and sequences in binary format
    """
    #first open the file outside "
    fin = open(fasta_name, 'rb')

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fin, lambda line: str(line, 'utf-8')[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = str(header.__next__(), 'utf-8')
        name=headerStr[1:].split(maxsplit=1)[0]
        seqlines= faiter.__next__()

        yield (name, seqlines)


def rename_fasta(fasta_in,fasta_out,map_names):

    assert type(map_names)==dict
    names= []

    with open(fasta_out,'wb') as fout:
        for name,seqlines in fasta_iter(fasta_in):

            if name in names:
                print(f"Seq with header {name} is duplicated")
            else:

                names.append(name)
                new_name=map_names.pop(name)

                fout.write(f">{new_name}\n".encode('utf-8'))
                fout.writelines(seqlines)








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
        print("ORF names are not unique, remove duplicates see:\n")

        print(orf2gene.loc[orf2gene.index.duplicated(keep=False)])

        orf2gene=orf2gene.loc[~orf2gene.index.duplicated()]

    orf2gene.to_csv(snakemake.output.cluster_attribution,sep='\t',header=True)

    # Rename representative sequence
    fasta_in=snakemake.input.faa
    fasta_out=snakemake.output.faa



    rename_fasta(fasta_in,fasta_out,map_names)
