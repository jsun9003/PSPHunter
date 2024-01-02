# -*- coding: utf-8 -*-
from sklearn.model_selection import LeaveOneOut
import numpy as np
from sklearn import svm
from scipy import stats
import random
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import roc_auc_score, f1_score, recall_score, precision_score, matthews_corrcoef, accuracy_score
import warnings
warnings.filterwarnings("ignore")

direct=sys.argv[1]
#matrix = np.loadtxt(direct + 'CVinput.txt',dtype=bytes).astype(str)
matrix = np.loadtxt(direct + 'CVinput.txt')
feature = np.array(matrix[:,1:])
tag = np.array(matrix[:,0])

re=10
kf = RepeatedKFold(n_splits=5, n_repeats=re,random_state=10086)
sum_auc = 0
'''
sum_recall = 0
sum_precision = 0
sum_f1 = 0
sum_acc = 0
sum_mcc = 0
'''
judge_final = np.zeros((101,6))
info = {}
for train, test in kf.split(feature):
	model_train = feature[train]
	model_test = feature[test]
	model_train_tag = tag[train]
	model_test_tag = tag[test]

	clf = RandomForestClassifier(n_estimators=500)
	clf.fit(model_train, model_train_tag) 

	predict_prob_y = clf.predict_proba(model_test)
	#print (test)
	for ind in range(0,len(test)):
		if test[ind] in info:
			sump=info[test[ind]]+np.array(predict_prob_y)[:,1][ind]
			info[test[ind]]=sump
		else:
			info[test[ind]]=np.array(predict_prob_y[:,1])[ind]
	test_auc = roc_auc_score(model_test_tag,predict_prob_y[:,1])
	sum_auc += test_auc
	#print test_auc

	recall = 0
	precision = 0
	f1 = 0
	acc = 0
	mcc = 0
	#prediction = clf.predict(model_test)
	judge = []
	
	for cutoff in range(0, 101):
		prediction = []
		judge_p = []
		for i in predict_prob_y[:, 1]:
			if i>cutoff*0.01:
				i=1
			else:
				i=0
			prediction.append(i)
		'print (len(prediction),len(model_test_tag))'
		'''
		sum_recall += recall_score(model_test_tag, prediction)
		sum_precision += precision_score(model_test_tag, prediction)
		sum_f1 += f1_score(model_test_tag, prediction)
		sum_acc += accuracy_score(model_test_tag, prediction)
		sum_mcc += matthews_corrcoef(model_test_tag, prediction)
		'''

		recall = recall_score(model_test_tag, prediction)
		precision = precision_score(model_test_tag, prediction)
		f1 = f1_score(model_test_tag, prediction)
		acc = accuracy_score(model_test_tag, prediction)
		mcc = matthews_corrcoef(model_test_tag, prediction)
		
		judge_p.append(0.01*cutoff)
		judge_p.append(recall)
		judge_p.append(precision)
		judge_p.append(f1)
		judge_p.append(acc)
		judge_p.append(mcc)
		judge.append(np.array(judge_p))
	#print (np.array(judge).shape,judge_final.shape)
	#print (judge[0],judge_final[0])
	judge_final += np.array(judge)
np.savetxt(direct + 'CVroc.txt', np.array(judge_final)/(5*re), fmt='%.3f',delimiter='\t')
with open(direct + 'CVroc.txt',"a") as f_out:
	f_out.write('%.3f\t' %(sum_auc/(5*re)))

f = open(direct + 'CVProb.txt','w')
#print (info)
for key in sorted(info.keys()):
	f.write('{:.3f}'.format((info[key])/re) + '\n')
f.close()
