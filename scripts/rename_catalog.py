import pandas as pd
import utils
import sys

import gzip as
sys.stdout= open(snakemake.log[0],"w")
sys.stderr= open(snakemake.log[0],"a")




def rename_fasta(fasta_in,fasta_out,new_names):

    old_names=[]
    n=0

    with open(fasta_in) as fin, gz.open(fasta_out,'wt') as fout:
        for line in fin:
            if line[0]=='>':

                old_name=line[1:].split(maxsplit=1)[0]

                assert (n==0) or (old_name!=old_names[-1]), f"Found duplicate representative name {old_name} in {fasta_in}"
                old_names.append(old_name)

                line=f">{new_names[n]} {old_name}\n"
                n+=1

            fout.write(line)
    return old_names

def parse_mmseqs_log_file(log_file,keyword="Number of clusters:"):

    with open(log_file) as f:

        result=None

        for line in f:
            if line.startswith(keyword):
                try:
                    result= int(line.replace(keyword,'').strip().rstrip())
                except ValueError as e:
                    raise Exception(f"Error parsing line:\n{line}\n") from e

        if result is None:
            raise Exception(f"Didn't found value in for keyword '{keyword}' in logfile {log_file}")
        else:
            return result


if __name__ == '__main__':



    Nrepresentatives = parse_mmseqs_log_file(snakemake.input.log)

    print(f"Number of representatives is {Nrepresentatives}")


    gene_names=utils.gen_names_for_range(Nrepresentatives,snakemake.params.prefix)

    original_names=  rename_fasta(snakemake.input.faa,snakemake.output.faa,gene_names)

    assert len(gene_names)==len(original_names), "Nuber of representatives should match the number of clusters found in the mmseqs log file"



    orf2gene = pd.Series(index=original_names,data=gene_names,name = 'Gene')
    orf2gene.index.name='ORF'



    orf2gene.to_csv(snakemake.output.name_mapping,sep='\t',header=True)

    del orf2gene

    # Rename representative sequence
