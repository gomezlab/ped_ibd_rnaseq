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


data = pd.read_csv('../data/clin_c.csv', index_col=0)

X = data.drop(['B2', 'B3', 'surg', 'rem', 'rsi', 'batch', 'mir'], axis=1)
y = data['surg']

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

keras_scores = pd.DataFrame(columns=['n_hidden', 'n_neurons', 'dropout', 'learning_rate', 'auroc', 'auprc'])

fold_cols = {}
count = 0
for n_hidden in [1, 2, 4]:
    for n_neurons in [100, 500, 1000]:
        for dropout in [0.2, 0.4, 0.8]:
            for learning_rate in [3e-3, 3e-4]:

                keras_y_true, keras_y_pred = [], []
        
                for train_index, val_index in cv.split(X):
            
                    X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                    y_train, y_val = y.iloc[train_index], y.iloc[val_index]

                    #feature selection with LASSO
                    lasso = LogisticRegression(penalty='l1', C=5, solver='liblinear', max_iter=9000)
                    lasso.fit(X_train, y_train)
                    model = SelectFromModel(lasso, prefit=True, max_features=20)
                    cols = model.get_support(indices=True)

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
        keras_scores.loc[count] = [n_hidden, n_neurons, dropout, learning_rate, auroc, auprc]
        keras_scores.to_csv('scores_nn_surg.csv')
        count += 1

keras_scores.to_csv('scores_nn_surg.csv')

xgb_scores = pd.DataFrame(columns=['n_estimators', 'max_depth', 'auroc', 'auprc'])


fold_cols = {}
count = 0

for n_estimators in [200, 500, 1000, 2000]:
    for max_depth in [4, 6, 12, 18]:
        
        xgb_y_true, xgb_y_pred = [], []

        for train_index, val_index in cv.split(X):
        
            X_train, X_val = X.iloc[train_index], X.iloc[val_index]
            y_train, y_val = y.iloc[train_index], y.iloc[val_index]

            #feature selection with LASSO
            lasso = LogisticRegression(penalty='l1', C=5, solver='liblinear', max_iter=9000)
            lasso.fit(X_train, y_train)
            model = SelectFromModel(lasso, prefit=True, max_features=20)
            cols = model.get_support(indices=True)            

            X_train = X_train.iloc[:, cols]
            X_val = X_val.iloc[:, cols]
            
            xgb = XGBClassifier(n_estimators=n_estimators, max_depth=max_depth, use_label_encoder=False, eval_metric='logloss', tree_method='gpu_hist', n_jobs=-1)
            xgb.fit(X_train, y_train)
            xgb_pred = xgb.predict_proba(X_val)
            xgb_y_pred.append(xgb_pred[:,1])
            xgb_y_true.append(y_val)
        
            del xgb

        #calculate auroc and auprc for this fold
        auroc = roc_auc_score(xgb_y_true, xgb_y_pred)
        auprc = average_precision_score(xgb_y_true, xgb_y_pred)
        xgb_scores.loc[count] = [n_estimators, max_depth, auroc, auprc]
        count += 1



xgb_scores.to_csv('scores_xgb_surg.csv')


rf_scores = pd.DataFrame(columns=['n_estimators', 'max_depth', 'auroc', 'auprc'])

fold_cols = {}
count = 0
for n_estimators in [200, 500, 1000, 2000]:
    for max_depth in [20, 60, 100, 200]:

        rf_y_true, rf_y_pred = [], []
        
        for train_index, val_index in cv.split(X):
            
                X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                y_train, y_val = y.iloc[train_index], y.iloc[val_index]

                #feature selection with LASSO
                lasso = LogisticRegression(penalty='l1', C=5, solver='liblinear', max_iter=9000)
                lasso.fit(X_train, y_train)
                model = SelectFromModel(lasso, prefit=True, max_features=20)
                cols = model.get_support(indices=True)

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
        rf_scores.loc[count] = [n_estimators, max_depth, auroc, auprc]
        count += 1

rf_scores.to_csv('scores_rf_surg.csv')

