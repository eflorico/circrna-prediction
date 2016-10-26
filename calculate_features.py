from parse import *
import random
import itertools
import time
import numpy as np
from pyfasta import Fasta
from intervaltree import IntervalTree
import sys

print("Loading files...")
t0 = time.time()

seqs = Fasta('data/hg19.fa')
crnas = read_bed('data/hsa_hg19_Rybak2015.bed')
exons = read_bed('data/all_exons.bed')
alus = read_bed('data/hg19_Alu.bed')
not_crnas = read_bed('tmp/negatives.bed')
not_crnas = not_crnas[:len(crnas)]

t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Calculating features...")
t0 = time.time()

# Concatenate data
data = crnas + not_crnas
labels = np.empty([ len(crnas) + len(not_crnas) ], np.int32)
labels[:len(crnas)] = 1
labels[len(crnas):] = 0

# Build ALU tree
alu_trees = {}
chromosomes = [ 'chr%d' % i for i in range(1, 23) ] + [ 'chrX', 'chrY' ]
for c in chromosomes:
	alu_trees[c][c.start:c.end] = True

# Build features
alu_flank_lengths = [ 50, 100, 200, 500 ]

K = 5
key_kmers = []
for i in range(1, K+1):
	key_kmers += [ "".join(c) for c in itertools.combinations_with_replacement('ACGT', i) ]

# Subset FIXME
subset = list(range(0, len(data)))
random.shuffle(subset)
subset = subset[0:len(data) // 10]

NUM_FEATURES = len(key_kmers) + 2 * len(alu_flank_lengths) + 1
features = np.empty([ len(subset), NUM_FEATURES ])
labels = labels[subset]

for i, j in enumerate(subset):
	gene = data[j]
	gene_data = get_gene_data(gene, seqs)

	# kmers
	kmer_features = [ 
		gene_data.count(kmer) / (gene.end - gene.start) 
		for kmer in key_kmers 
	]

	# ALU counts
	alu_counts = [
		len(alu_trees[gene.chr_n][gene.start - length:gene.start])
		for length in alu_flank_lengths
	] + [
		len(alu_trees[gene.chr_n][gene.end:gene.end + length])
		for length in alu_flank_lengths
	] + [
		alu_trees[gene.chr_n][gene.start:gene.end]
	]

	# Concatenate features
	features[i] = kmer_features + alu_counts

	if i % 1000 == 0:
		print('.', end='')
		sys.stdout.flush()

np.save('tmp/features.npy', features)
np.save('tmp/labels.npy', labels)

t1 = time.time()
print("%.2fs" % (t1 - t0))