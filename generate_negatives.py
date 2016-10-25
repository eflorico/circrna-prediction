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

t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Loading cRNA data...")
t0 = time.time()

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

bedFile = open('tmp/negatives.bed', 'w')

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

	# Write to file
	print("%s\t%9d\t%9d\t%s" % (chr_n, start, end, strand), file=bedFile)
bedFile.close()

t1 = time.time()
print("%.2fs" % (t1 - t0))