import re
from intervaltree import Interval, IntervalTree
from collections import namedtuple

Gene = namedtuple('Gene', ['chr_n', 'start', 'end', 'gene_name', 'strand'])

def read_fasta(file):
    from Bio import SeqIO

    fasta = SeqIO.parse(file, "fasta")

    # Only use chrN, chrX, chrY
    # chromosome_pattern = re.compile('^chr(\d+|X|Y)$')

    seqs = {}

    for record in fasta:
        #if chromosome_pattern.match(record.id) != None:
        # Make all upper case
        seqs[record.id] = record.seq.upper()
        print("%s sequence loaded" % record.id)
        break # FIXME

    fasta.close()

    return seqs

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

        # FIXME
        if chr_n != 'chr1':
            continue

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

def build_interval_trees(genes, seqs):
    trees = {}
    for seq_name, _ in seqs.items():
        trees[seq_name] = IntervalTree()
        
    for gene in genes:
        trees[gene.chr_n][gene.start:gene.end] = gene

    return trees