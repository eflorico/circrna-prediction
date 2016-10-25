from intervaltree import IntervalTree
import numpy as np
from sklearn import svm
from sklearn.model_selection import KFold
from parse import *
import random
import itertools
import time
from sklearn.metrics import roc_auc_score, precision_recall_curve, f1_score

print("Loading files...")
t0 = time.time()

seqs = read_fasta('data/hg19.fa')
crnas = read_bed('data/hsa_hg19_Rybak2015.bed')
exons = read_bed('data/all_exons.bed')

t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Loading cRNA data...")
t0 = time.time()

#crna_data = [ get_gene_data(crna, seqs) for crna in crnas if crna.chr_n == 'chr1' ]

t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Building trees...")
t0 = time.time()

# Build negative examples
exon_trees = build_interval_trees(exons, seqs)
free_trees = build_interval_trees([], seqs)
for chr_n, tree in free_trees.items():
	tree[exon_trees[chr_n].begin():exon_trees[chr_n].end()] = True

exon_trees['chr1'].remove_overlap(179079416, 179091002)

for crna in crnas:
	exon_trees[crna.chr_n].remove_overlap(crna.start, crna.end)
	free_trees[crna.chr_n].chop(crna.start, crna.end)

t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Generating negative examples...")
t0 = time.time()

not_crnas = []
for crna in crnas:
	# Choose length
	length = random.randint(100, 10000)

	# Use same chromosome distribution FIXME?
	start, end = crna.start, crna.end
	chr_n = crna.chr_n

	# Get trees for chromosome sequence
	ft = free_trees[chr_n]
	et = exon_trees[chr_n]

	# Find spot with exactly that length
	free_intvs = list(ft.items())

	while True:
		intv = random.choice(free_intvs)

		if intv.length() < length: 
			continue

		start = intv.begin + random.randint(0, intv.length() - length)
		end = start + length

		if et.overlaps(start) and et.overlaps(end):
			break

	# Choose strand
	strand = random.choice('+-')
	gene = Gene(chr_n, start, end, "neg sample", strand)
	#not_crnas.append(get_gene_data(gene, seqs))
	not_crnas.append(gene)

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
print(key_kmers)

kmer_features = np.empty([ len(data), len(key_kmers) ])
for i, gene in enumerate(data):
	gene_data = get_gene_data(gene, seqs)
	kmer_features[i] = [ gene_data.count(kmer) / (gene.end - gene.start) 
		for kmer in key_kmers ]

# SVM
t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Testing SVM...")
t0 = time.time()

clf = svm.SVC()
kf = KFold(n_splits=10, shuffle=True)

print("ROC AUC\tF1")

for train, test in kf.split(kmer_features):
	clf.fit(kmer_features[train], labels[train])

	y_true = labels[test]
	y_pred = clf.predict(kmer_features[test])

	roc_auc = roc_auc_score(y_true, y_pred)
	f1 = f1_score(y_true, y_pred)

	print("%.4f\t%.4f", roc_auc, f1)

t1 = time.time()
print("%.2fs" % (t1 - t0))