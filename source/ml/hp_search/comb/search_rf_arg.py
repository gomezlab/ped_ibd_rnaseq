import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--outcome', type=str, default='B2')
parser.add_argument('--rsi', type=str, default='no_rsi')

args = parser.parse_args()
outcome = args.outcome
rsi = args.rsi

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve, accuracy_score, auc, precision_recall_curve, average_precision_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_score
from sklearn.metrics import recall_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from xgboost.sklearn import XGBClassifier
from matplotlib import pyplot as plt
from math import sqrt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from scipy import stats as st
from random import randrange
from sklearn.model_selection import LeaveOneOut
from tensorflow import keras
from sklearn.feature_selection import SelectFromModel


data = pd.read_csv('../../data/comb_vst_c.csv', index_col=0)

clin_cols = ['sex', 'dx_age', 'initial_b', 'location', 'perianal_disease']

X = data.drop(['B2', 'B3', 'surg', 'rem', 'rsi'], axis=1)
if rsi == 'rsi':
    X['rsi'] = data['rsi']
    clin_cols.append('rsi')
y = data[outcome]

cv = LeaveOneOut()

rf_scores = pd.DataFrame(columns=['feature_count', 'n_estimators', 'max_depth', 'auroc', 'auprc'])

fold_cols = {}
count = 0
for feat_count in [0.1, 1, 5, 10]:
    for n_estimators in [200, 500, 1000, 2000]:
        for max_depth in [20, 60, 100, 200]:

            rf_y_true, rf_y_pred = [], []
            
            for train_index, val_index in cv.split(X):
                
                    X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                    y_train, y_val = y.iloc[train_index], y.iloc[val_index]

                    #feature selection with LASSO
                    lasso = LogisticRegression(penalty='l1', C=feat_count, solver='liblinear', random_state=0, max_iter=9000)
                    lasso.fit(X_train, y_train)
                    model = SelectFromModel(lasso, prefit=True)
                    cols = model.get_support(indices=True)
                    for clin_var in clin_cols:
                        clin_idx = X_train.columns.get_loc(clin_var)
                        if clin_idx not in cols:
                            cols = np.append(cols, clin_idx)
                        else:
                            continue

                    X_train = X_train.iloc[:, cols]
                    X_val = X_val.iloc[:, cols]

                    rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, n_jobs=-1)
                    rf.fit(X_train, y_train)
                    rf_pred = rf.predict_proba(X_val)
                    rf_y_pred.append(rf_pred[:,1])
                    rf_y_true.append(y_val)

                    del rf

            #calculate auroc and auprc for this fold
            auroc = roc_auc_score(rf_y_true, rf_y_pred)
            auprc = average_precision_score(rf_y_true, rf_y_pred)
            rf_scores.loc[count] = [feat_count, n_estimators, max_depth, auroc, auprc]
            count += 1

rf_scores.sort_values(by='auroc', ascending=False, inplace=True)
rf_scores.to_csv('results/{}/rf_{}.csv'.format(rsi, outcome))

