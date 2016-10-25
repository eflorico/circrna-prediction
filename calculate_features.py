from intervaltree import IntervalTree
from sklearn import svm
from sklearn.model_selection import KFold
from parse import *
import random
import itertools
import time
import numpy as np
from sklearn.metrics import roc_auc_score, precision_recall_curve, f1_score

print("Loading files...")
t0 = time.time()

seqs = read_fasta('data/hg19.fa')
crnas = read_bed('data/hsa_hg19_Rybak2015.bed')
exons = read_bed('data/all_exons.bed')
not_crnas = read_bed('tmp/negatives.bed')

t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Calculating features...")
t0 = time.time()

# Concatenate data
data = crnas + not_crnas
labels = np.empty([ len(crnas) + len(not_crnas) ], np.int32)
labels[:len(crnas)] = 1
labels[len(crnas):] = 0

# Build features
k = 3

key_kmers = []
for i in range(1, k+1):
	key_kmers += [ "".join(c) for c in itertools.combinations_with_replacement('ACGT', i) ]

kmer_features = np.empty([ len(data), len(key_kmers) ])
for i, gene in enumerate(data):
	gene_data = get_gene_data(gene, seqs)
	kmer_features[i] = [ gene_data.count(kmer) / (gene.end - gene.start) 
		for kmer in key_kmers ]
np.save('tmp/features.npy', kmer_features)
np.save('tmp/labels.npy', labels)

t1 = time.time()
print("%.2fs" % (t1 - t0))