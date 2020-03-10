import os

def _make_test_fasta(test_file='test_ABC.fasta'):
    with open(test_file,'w') as f:
        for Number,Letter in enumerate("ATCG"):
            f.write(f'>contig_{Number+1} description\n{Letter}\n')

def str2multiline(long_string,line_length=80):
    return [long_string[i:i+line_length]   for i in range(0,len(long_string),line_length)]

def count_Nseq(fasta_file):
    """
    Counts number of sequences in a fasta file.
    >>> fasta_file='test_ABC.fasta'
    >>> _make_test_fasta(fasta_file) # makes fasta with a seq for each nucleotide
    >>> count_Nseq(fasta_file)
    4
    >>> os.remove(fasta_file)

    """
    i=0
    with open(fasta_file) as f:
        for line in f:
            if line[0]=='>':
                i+=1
    return i

def split(fasta_file,SubsetSize,out_dir):
    """
    Splits a fasta in subsets of size max maxSubsetSize.
    >>> fasta_file='test_ABC.fasta'
    >>> out_dir = 'test_outdit_doctest'
    >>> _make_test_fasta(fasta_file) # makes fasta with a seq for each nucleotide
    >>> split(fasta_file,3,out_dir)
    >>> len(os.listdir(out_dir))
    2
    >>> count_Nseq('test_outdit_doctest/subset1.fasta')
    3
    >>> count_Nseq('test_outdit_doctest/subset2.fasta')
    1
    >>> split(fasta_file,3,out_dir)
    Traceback (most recent call last):
        ...
    FileExistsError: [Errno 17] File exists: 'test_outdit_doctest'
    >>> import shutil; shutil.rmtree(out_dir)
    >>> os.remove(fasta_file)
    """
    extension= os.path.splitext(fasta_file)[-1]
    os.makedirs(out_dir)


    i,subset_n=0,0
    fout= None

    with open(fasta_file,'r') as fin:
        for line in fin:

            if line[0]=='>':

                if (i % SubsetSize) == 0:
                    subset_n+=1
                    if fout is not None:
                        fout.close()
                    fout = open(f"{out_dir}/subset{subset_n}{extension}",'w')

                i+=1
            fout.write(line)

def header2origin(fasta_file,out,simplify_header=True):
    """
    Annotates a fasta file to it's filename:
    genome.fasta:
    >contig1 description
    ACTAC
    >contig2 description
    ACTAC
    ...

    becomes:
    contig1 genome
    contig2 genome

    input is a fasta filename:
    out is a filename or a stream

    """

    if type(out) == str:
        out_stream=open(out,'w')
    else:
        out_stream= out

    name= os.path.splitext(os.path.split(fasta_file)[-1])[0]

    # write names of contigs in mapping file
    with open(fasta_file) as f :
        for line in f:
            if line[0]==">":
                header=line[1:].strip()
                if simplify_header:
                    header=header.split()[0]
                out_stream.write(f"{header}\t{name}\n")
    out_stream.flush()

from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of name as string and sequences
    """
    #first open the file outside "
    fin = open(fasta_name, 'r')

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fin, lambda line: line[0] == ">"))

    for header in faiter:
        headerStr = header.__next__()
        name=headerStr[1:].split(maxsplit=1)[0]
        seqlines= faiter.__next__()



        yield (name,headerStr ,seqlines)


if __name__ == "__main__":
    import doctest, shutil
    doctest.testmod()
