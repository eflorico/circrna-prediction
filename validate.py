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
kmer_features = np.load('tmp/features.npy')
labels = np.load('tmp/labels.npy')

# SVM
t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Testing SVM...")
t0 = time.time()

clf = svm.SVC()
kf = KFold(n_splits=10, shuffle=True)

print("%6s %6s" % ("ROC", "F1"))

for train, test in kf.split(kmer_features):
	clf.fit(kmer_features[train], labels[train])

	y_true = labels[test]
	y_pred = clf.predict(kmer_features[test])

	roc_auc = roc_auc_score(y_true, y_pred)
	f1 = f1_score(y_true, y_pred)

	print("%.4f\t%.4f" % (roc_auc, f1))

t1 = time.time()
print("%.2fs" % (t1 - t0))

