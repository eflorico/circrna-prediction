import re
import time
from collections import namedtuple

Gene = namedtuple('Gene', [
    'chr_n', 
    'start', 
    'end', 
    'gene_name', 
    'strand', 
    'blocks'
])

def tick(label=None):
    t1 = time.time()

    if tick.t0 != None:
        print("%.4fs" % (t1 - tick.t0))

    if label != None:
        print(label)

    tick.t0 = t1

tick.t0 = None

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
    chr_filter = re.compile('^chr(\d+|X|Y)\s')
    genes = []

    for line in bedFile:
        # Filter chromosomes
        if chr_filter.match(line) == None:
            continue

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

        # Very long BED
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

        gene = Gene(chr_n, int(start), int(end), gene_name, strand, blocks)
        genes.append(gene)

    bedFile.close()
    return genes

def group_by_chromosome(genes):
    groups = {}

    for gene in genes:
        groups.setdefault(gene.chr_n, []).append(gene)

    return groups

def get_gene_data(gene, seqs):
    data = seqs[gene.chr_n][gene.start:gene.end].upper()

    if gene.strand == '-':
        data = reverse_complement(data)

    return data