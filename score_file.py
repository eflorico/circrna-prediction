from sklearn.metrics import roc_auc_score, precision_recall_curve, f1_score, average_precision_score

# Score prediction from lsgkm

results = open('lsgkm/bin/crna.cvpred.txt', 'r')
y_true, y_pred = [], []

for line in results:
	values = line.split("\t")
	y_true.append(float(values[2]))
	y_pred.append(float(values[1]))

results.close()

roc_auc = roc_auc_score(y_true, y_pred)
aupr = average_precision_score(y_true, y_pred)

y_pred_bin = np.zeros([ len(y_pred) ])
y_pred_bin[:] = y_pred
y_pred_bin = y_pred_bin > 0
f1 = f1_score(y_true, y_pred_bin)
f1=0

print("%.4f\t%.4f\t%.4f\t" % (roc_auc, aupr, f1))
