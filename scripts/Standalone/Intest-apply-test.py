# -*- coding: utf-8 -*-
import numpy as np
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, f1_score, recall_score, precision_score, matthews_corrcoef, accuracy_score
import joblib
import warnings
warnings.filterwarnings("ignore")

direct1=sys.argv[1]
direct2=sys.argv[2]
re=sys.argv[3]
#matrix = np.loadtxt(direct + 'CVinput.txt',dtype=bytes).astype(str)

        
clf = joblib.load(direct1 + "train_model.m")
matrixI = np.loadtxt(direct2 + 'Intestinput.txt')
matrixI = matrixI.reshape(matrixI.shape[0], -1)

featureI = np.array(matrixI[:,1:])
tagI = np.array(matrixI[:,0])

model_test = featureI
predict_prob_y = clf.predict_proba(model_test)

f = open(direct2 + 'InProbphase' + re + '.txt','w')
for i in range(len(predict_prob_y)):
	f.write('{:.3f}\n'.format(predict_prob_y[:,1][i]))
f.close()
