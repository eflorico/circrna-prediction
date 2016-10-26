import re
from intervaltree import Interval, IntervalTree
from collections import namedtuple

Gene = namedtuple('Gene', ['chr_n', 'start', 'end', 'gene_name', 'strand'])

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

        # Short BED
        if vals[3] in ['+','-']:
            strand = vals[3]
            gene_name = ""
        # Long BED
        else:
            strand = vals[5]
            gene_name = vals[3]

        start, end = int(start), int(end)

        # Omit weird data
        if end <= start: continue

        genes.append(Gene(chr_n, int(start), int(end), gene_name, strand))

    bedFile.close()
    return genes

def get_gene_data(gene, seqs):
    data = seqs[gene.chr_n][gene.start:gene.end]

    if gene.strand == '-':
        data = reverse_complement(data)

    return data

ch_filter = re.compile('chr(\d+|X|Y)$')

def build_interval_trees(genes, seqs):
    trees = {}
    for seq_name, _ in seqs.items():
        if ch_filter.match(seq_name) != None:
            trees[seq_name] = IntervalTree()
        
    for gene in genes:
        if ch_filter.match(gene.chr_n) != None:
            trees[gene.chr_n][gene.start:gene.end] = gene

    return trees