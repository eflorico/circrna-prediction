import time
import numpy as np
from sklearn import svm
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, precision_recall_curve, f1_score

features = np.load('tmp/features.npy')
labels = np.load('tmp/labels.npy')

# k Fold
def test(model, params):
	print(', '.join('%s: %s' % (k,v) for k, v in params.items()))

	kf = KFold(n_splits=10, shuffle=True)

	scores = []
	for train, test in kf.split(features):
		try:
			model.fit(features[train], labels[train])
		except KeyboardInterrupt:
			print("Aborted")
			return
		finally:
			y_true = labels[test]
			y_pred = clf.predict(features[test])

			roc_auc = roc_auc_score(y_true, y_pred)
			f1 = f1_score(y_true, y_pred)
			scores.append([ roc_auc, f1 ])

			print("%.4f\t%.4f" % (roc_auc, f1))

	if len(scores) > 0:
		print("Mean scores:")
		n = len(scores)
		m = len(scores[0])
		means = [ 
			sum(s[i] for s in scores) / n
			for i in range(0, m) 
		]
		print(("%.4f\t" * m) % tuple(means))

# Random Forest
print("Testing RF...")
t0 = time.time()

params = { 'n_estimators': 120 }
clf = RandomForestClassifier(**params)
test(clf, params)

# SVM
t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Testing SVM...")
t0 = time.time()

params = { 'kernel': 'rbf' }
clf = svm.SVC(**params)
test(clf, params)

t1 = time.time()
print("%.2fs" % (t1 - t0))