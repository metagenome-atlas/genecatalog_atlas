import pandas as pd


def get_names(fasta_in):

    names=[]

    with open(fasta_in) as fin:
        for line in fin:
            if line[0]=='>':
                name=line[1:].split(maxsplit=1)[0]
                names.append(name)
    return names

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

fasta_in=snakemake.input.faa
fasta_out=snakemake.output.faa

mapping_file= snakemake.input.mapping

representatives= get_names(fasta_in)

assert len(representatives)==len(set(representatives))

Mapping=pd.read_csv(mapping_file,index_col=0,sep='\t',squeeze=True)
map_names= Mapping.loc[representatives].to_dict()
del Mapping

rename_fasta(fasta_in,fasta_out,map_names)
