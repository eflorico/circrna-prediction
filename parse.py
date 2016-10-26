import re
from collections import namedtuple

Gene = namedtuple('Gene', [
    'chr_n', 
    'start', 
    'end', 
    'gene_name', 
    'strand', 
    'blocks'
])

def complement(seq):
    complement = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' }
    complseq = [complement[base] for base in seq]
    return complseq

def reverse_complement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

def read_bed(file):
    bedFile = open(file, 'r')
    genes = []

    for line in bedFile:
        vals = line.split()
        chr_n, start, end = vals[0:3]
        blocks = []

        # Short BED
        if vals[3] in ['+','-']:
            strand = vals[3]
            gene_name = ""
            
        # Long BED
        else:
            strand = vals[5]
            gene_name = vals[3]

        if len(vals) > 11:
            block_starts = [ int(x) for x in vals[11].split(',') ]
            block_lens   = [ int(x) for x in vals[10].split(',') ]
            blocks = [ 
                (block_starts[i], block_lens[i]) 
                for i in range(0, len(block_starts))
            ]

        start, end = int(start), int(end)

        # Omit weird data
        if end <= start: continue

        genes.append(Gene(chr_n, int(start), int(end), gene_name, strand, blocks))

    bedFile.close()
    return genes

def get_gene_data(gene, seqs):
    data = seqs[gene.chr_n][gene.start:gene.end].upper()

    if gene.strand == '-':
        data = reverse_complement(data)

    return data