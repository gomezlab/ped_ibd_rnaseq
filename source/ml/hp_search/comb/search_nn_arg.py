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

def build_model(n_hidden=4, n_neurons=50, dropout=0.2, learning_rate=3e-3):
    model = keras.models.Sequential()
    model.add(keras.layers.Flatten(input_shape=input_shape))
    model.add(keras.layers.BatchNormalization())
    for layer in range(n_hidden):
        model.add(keras.layers.Dense(n_neurons))
        model.add(keras.layers.BatchNormalization())
        model.add(keras.layers.Dropout(dropout))
        model.add(keras.layers.Activation("relu"))
    model.add(keras.layers.Dense(1, activation="sigmoid"))

    optimizer = keras.optimizers.Adam(learning_rate=learning_rate)
    
    model.compile(loss=keras.losses.BinaryCrossentropy(), metrics=['AUC'], optimizer=optimizer)
    return model

keras_scores = pd.DataFrame(columns=['feature_count', 'n_hidden', 'n_neurons', 'dropout', 'learning_rate', 'auroc', 'auprc'])

fold_cols = {}
count = 0

for feature_count in [0.5, 1, 5, 10]:
    for n_hidden in [1, 2, 4]:
        for n_neurons in [100, 500, 1000]:
            for dropout in [0.2, 0.4, 0.8]:
                for learning_rate in [3e-3, 3e-4]:

                    keras_y_true, keras_y_pred = [], []
            
                    for train_index, val_index in cv.split(X):
                
                        X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                        y_train, y_val = y.iloc[train_index], y.iloc[val_index]

                        #feature selection with LASSO
                        lasso = LogisticRegression(penalty='l1', C=feature_count, solver='liblinear', random_state=0, max_iter=9000)
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

                        input_shape = X_train.shape[1:]

                        keras_obj = lambda: build_model(n_hidden=n_hidden, n_neurons=n_neurons, dropout=dropout, learning_rate=learning_rate)

                        keras_clf = keras.wrappers.scikit_learn.KerasClassifier(build_fn=keras_obj)
                        keras_clf.fit(X_train, y_train, epochs=150, batch_size=32, validation_data=(X_val, y_val), verbose=0)
                        keras_pred = keras_clf.predict_proba(X_val)
                        keras_y_pred.append(keras_pred[:,1])
                        keras_y_true.append(y_val)

                        del keras_clf

                    #calculate auroc and auprc for this fold
                    auroc = roc_auc_score(keras_y_true, keras_y_pred)
                    auprc = average_precision_score(keras_y_true, keras_y_pred)
                    keras_scores.loc[count] = [feature_count, n_hidden, n_neurons, dropout, learning_rate, auroc, auprc]
                    count += 1
                    del keras_y_true, keras_y_pred
                    
keras_scores.sort_values(by='auroc', ascending=False, inplace=True)
keras_scores.to_csv('results/{}/nn_{}.csv'.format(rsi, outcome))



