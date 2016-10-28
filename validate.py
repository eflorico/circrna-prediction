import time
import itertools
import random
import numpy as np
from util import tick
from sklearn import svm
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import roc_auc_score, precision_recall_curve, f1_score, average_precision_score

features = np.load('tmp/features.npy')
labels = np.load('tmp/labels.npy')

# Zero-center
features = (features - features.mean(axis=0))

# Subset for grid search
random.seed(0)
subset = list(range(0, len(features)))
subset = random.sample(subset, 1000)

# k Fold
def test(model):
	print(', '.join('%s: %s' % (k,v) for k, v in model.get_params().items()))
	print("%6s\t%6s\t%6s" % ('AUC', 'AUPR', 'F1'))

	kf = KFold(n_splits=10, shuffle=True)

	scores = []
	for train, test in kf.split(features):
		try:
			model.fit(features[train], labels[train])
		
			y_true = labels[test]
			y_pred = np.array(clf.predict_proba(features[test]))[:, 1]

			roc_auc = roc_auc_score(y_true, y_pred)
			aupr = average_precision_score(y_true, y_pred)
			f1 = f1_score(y_true, y_pred >= .5)
			scores.append([ roc_auc, aupr, f1 ])

			print("%.4f\t%.4f\t%.4f" % (roc_auc, aupr, f1))
		except KeyboardInterrupt:
			print("Aborted")
			return

	if len(scores) > 0:
		print("Mean scores:")
		n = len(scores)
		m = len(scores[0])
		means = [ 
			sum(s[i] for s in scores) / n
			for i in range(0, m) 
		]
		print(("%.4f\t" * m) % tuple(means))

tick("Random Forest...")
params = { 'n_estimators': 120, 'n_jobs': -1 }
clf = RandomForestClassifier(**params)
test(clf)
exit()

# tick("Feature ranking...")

# # Fields
# alu_flank_lengths = [ 50, 100, 200, 500 ]
# K = 7

# # kmer field names
# feature_names = []
# for i in range(1, K+1):
# 	feature_names += [ "".join(c) for c in itertools.combinations_with_replacement('ACGT', i) ]

# for length in alu_flank_lengths:
# 	feature_names.extend([
# 		'-%d/0 left ALUs' % length,
# 		'-%d/+%d left ALUs' % (length, length),
# 		'0/+%d left ALUs' % length,
# 		'-%d/0 right ALUs' % length,
# 		'-%d/+%d right ALUs' % (length, length),
# 		'0/+%d right ALUs' % length
# 	])

# clf.fit(features, labels)
# ranking = sorted(zip(clf.feature_importances_, feature_names), reverse=True)
# for rank, name in ranking:
# 	print("%.4f %s" % (rank, name))


tick("Lin Reg...")
params = { 'n_jobs': -1, 'C': 10000 }
param_grid = { 'n_jobs': [-1], 'C': [ 10 ** i for i in range(4, 10)] }
clf = GridSearchCV(LogisticRegression(), param_grid)
clf.fit(features[subset], labels[subset])
print(clf.best_estimator_)
test(clf.best_estimator_)

# tick("kNN")
# clf = KNeighborsClassifier(n_jobs=-1)
# test(clf)

tick("SVM grid search...")
params = { 'C': 2 ** 6, 'gamma': 2 ** -9, 'kernel': 'rbf' }
param_grid = [
#	{ 'kernel':['rbf']}
   { 'C': [10**i for i in range(3,8)], 'kernel': ['linear']},
   { 'C': [2 ** 6], 'gamma': [2 ** -9], 'kernel': ['rbf'] }
#  { 'C': [ 10 ** i for i in range(4, 9) ], 
#  	'gamma': [ 2 ** i for i in range(-10, -4) ], 
#  	'kernel': ['rbf']
#  },
# {'C': [1000, 10000], 'degree': [3, 4], 'kernel': ['poly']},
 ]

#train, test = KFold(n_splits=100, shuffle=True)

clf = GridSearchCV(svm.SVC(), param_grid)
clf.fit(features[subset], labels[subset])
print(clf.best_estimator_)

tick("SVM...")
params = { 'kernel': 'rbf' }
test(clf.best_estimator_)

tick()