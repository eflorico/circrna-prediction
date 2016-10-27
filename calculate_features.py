from util import *
import random
import itertools
import time
import numpy as np
from pyfasta import Fasta
from intervaltree import IntervalTree
import sys

tick("Loading files...")

seqs = Fasta('data/hg19.fa')
crnas = read_bed('data/hsa_hg19_Rybak2015.bed')
#exons = read_bed('data/all_exons.bed')
alus = group_by_chromosome(read_bed('data/hg19_Alu.bed'))
not_crnas = read_bed('tmp/negatives.bed')

# Don't use more negatives then necessary
not_crnas = not_crnas[:len(crnas)]

tick("Preparing features...")

# Concatenate data
data = crnas + not_crnas
labels = np.empty([ len(data) ], np.int32)
labels[:len(crnas)] = 1
labels[len(crnas):] = 0

# Build ALU tree
alu_trees = {}
for chr_n, genes in alus.items():
	alu_trees[chr_n] = IntervalTree()

	for alu in genes:
		alu_trees[chr_n][alu.start:alu.end] = True

tick("Building features...")

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
		len(alu_trees[gene.chr_n][gene.start:gene.end])
	]

	# Concatenate features
	features[i] = kmer_features + alu_counts

	if i % 1000 == 0:
		print('.', end='')
		sys.stdout.flush()

np.save('tmp/features.npy', features)
np.save('tmp/labels.npy', labels)

tick()