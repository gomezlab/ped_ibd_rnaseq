#%%
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

#%%
for outcome in ['B2', 'B3', 'rem', 'surg']:
    #%%
    data = pd.read_csv('../data/clin_c.csv', index_col=0)
    rf_params = pd.read_csv('scores_rf_{}.csv'.format(outcome), index_col=0)

    rf_params.head()
    #%%
    rf_params.reset_index(inplace=True, drop=True)
    
    rf_n_estimators = (rf_params['n_estimators'][0]).astype(int)
    rf_max_depth = (rf_params['max_depth'][0]).astype(int)

    #%%


    X = data.drop(['B2', 'B3', 'surg', 'rem', 'batch', 'mir'], axis=1)
    y = data[outcome]

    cv = LeaveOneOut()
    rf_y_true, rf_y_pred = [], []

    
    for n_estimators in [rf_n_estimators]:
        for max_depth in [rf_max_depth]:


            
            for train_index, val_index in cv.split(X):
                
                    X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                    y_train, y_val = y.iloc[train_index], y.iloc[val_index]

                    rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, n_jobs=-1)
                    rf.fit(X_train, y_train)
                    rf_pred = rf.predict_proba(X_val)
                    rf_y_pred.append(rf_pred[:,1])
                    rf_y_true.append(y_val)

                    del rf


    data = pd.read_csv('../data/clin_c.csv', index_col=0)

    xgb_params = pd.read_csv('scores_xgb_{}.csv'.format(outcome), index_col=0)
    xgb_params.reset_index(inplace=True, drop=True)
    
    xgb_n_estimators = (xgb_params['n_estimators'][0]).astype(int)
    xgb_max_depth = (xgb_params['max_depth'][0]).astype(int)

    xgb_y_true, xgb_y_pred = [], []


    X = data.drop(['B2', 'B3', 'surg', 'rem', 'batch', 'mir'], axis=1)
    y = data[outcome]
    
    for n_estimators in [xgb_n_estimators]:
        for max_depth in [xgb_max_depth]:        

            for train_index, val_index in cv.split(X):
            
                X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                y_train, y_val = y.iloc[train_index], y.iloc[val_index]             
         
                xgb = XGBClassifier(n_estimators=n_estimators, max_depth=max_depth, use_label_encoder=False, eval_metric='logloss', tree_method='gpu_hist', n_jobs=-1)
                xgb.fit(X_train, y_train)
                xgb_pred = xgb.predict_proba(X_val)
                xgb_y_pred.append(xgb_pred[:,1])
                xgb_y_true.append(y_val)
            
                del xgb


    data = pd.read_csv('../data/clin_c.csv', index_col=0)
    nn_params = pd.read_csv('scores_nn_{}.csv'.format(outcome), index_col=0)
    nn_params.reset_index(inplace=True, drop=True)
    
    nn_n_hidden = (nn_params['n_hidden'][0]).astype(int)
    nn_n_neurons = (nn_params['n_neurons'][0]).astype(int)
    nn_dropout = nn_params['dropout'][0]
    nn_learning_rate = nn_params['learning_rate'][0]

    X = data.drop(['B2', 'B3', 'surg', 'rem', 'batch', 'mir'], axis=1)
    y = data[outcome]

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

    cv = LeaveOneOut()

    for n_hidden in [nn_n_hidden]:
        for n_neurons in [nn_n_neurons]:
            for dropout in [nn_dropout]:
                for learning_rate in [nn_learning_rate]:

                    keras_y_true, keras_y_pred = [], []
            
                    for train_index, val_index in cv.split(X):
                
                        X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                        y_train, y_val = y.iloc[train_index], y.iloc[val_index]


                        input_shape = X_train.shape[1:]

                        keras_obj = lambda: build_model(n_hidden=n_hidden, n_neurons=n_neurons, dropout=dropout, learning_rate=learning_rate)

                        keras_clf = keras.wrappers.scikit_learn.KerasClassifier(build_fn=keras_obj)
                        keras_clf.fit(X_train, y_train, epochs=150, batch_size=32, validation_data=(X_val, y_val), verbose=0)
                        keras_pred = keras_clf.predict_proba(X_val)
                        keras_y_pred.append(keras_pred[:,1])
                        keras_y_true.append(y_val)

                        del keras_clf

    cv = LeaveOneOut()
    lr_y_true, lr_y_pred = [], []

            
    for train_index, val_index in cv.split(X):

        X_train, X_val = X.iloc[train_index], X.iloc[val_index]
        y_train, y_val = y.iloc[train_index], y.iloc[val_index]

        lr = LogisticRegression()
        lr.fit(X_train, y_train)
        lr_pred = lr.predict_proba(X_val)
        lr_y_pred.append(lr_pred[:,1])
        lr_y_true.append(y_val)

        del lr

    def auroc_ci(y_true, y_pred):
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        roc_auc = auc(fpr, tpr)
        mean = roc_auc
        std = sqrt(roc_auc * (1.0 - roc_auc) / len(y_true))
        low  = mean - std
        high = mean + std
        return low, mean, high
    rf_confidence = auroc_ci(rf_y_true, rf_y_pred)
    xgb_confidence = auroc_ci(xgb_y_true, xgb_y_pred)
    keras_confidence = auroc_ci(keras_y_true, keras_y_pred)
    lr_confidence = auroc_ci(lr_y_true, lr_y_pred)

    #create labels for roc curves
    rf_label = 'RF: AUROC ' + str(round(rf_confidence[1], 3)) + ' (95% CI ' + str(round(rf_confidence[0], 3)) + ' - ' + str(round(rf_confidence[2], 3)) + ')'
    xgb_label = 'XGB: AUROC ' + str(round(xgb_confidence[1], 3)) + ' (95% CI ' + str(round(xgb_confidence[0], 3)) + ' - ' + str(round(xgb_confidence[2], 3)) + ')'
    keras_label = 'NN: AUROC ' + str(round(keras_confidence[1], 3)) + ' (95% CI ' + str(round(keras_confidence[0], 3)) + ' - ' + str(round(keras_confidence[2], 3)) + ')'
    lr_label = 'LASSO: AUROC ' + str(round(lr_confidence[1], 3)) + ' (95% CI ' + str(round(lr_confidence[0], 3)) + ' - ' + str(round(lr_confidence[2], 3)) + ')'

    #calculate tpr and fpr for each model
    rf_fpr, rf_tpr, _ = roc_curve(rf_y_true, rf_y_pred)
    xgb_fpr, xgb_tpr, _ = roc_curve(xgb_y_true, xgb_y_pred)
    keras_fpr, keras_tpr, _ = roc_curve(keras_y_true, keras_y_pred)
    lr_fpr, lr_tpr, _ = roc_curve(lr_y_true, lr_y_pred)

    import matplotlib
    matplotlib.rcParams.update({'font.size': 16})
    #plot the ROC curves for each model
    plt.figure(figsize=(10,10))
    plt.plot(lr_fpr, lr_tpr, color='cadetblue', label=lr_label)
    plt.plot(rf_fpr, rf_tpr, color='deepskyblue', label=rf_label)
    plt.plot(xgb_fpr, xgb_tpr, color='steelblue', label=xgb_label)
    plt.plot(keras_fpr, keras_tpr, color='dodgerblue', label=keras_label)
    plt.plot([0, 1], [0, 1], color='black', linestyle='--')
    plt.legend(loc="lower right")
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.xlim([0.0, 1.0])
    plt.savefig('results/roc_clin_{}_norect.png'.format(outcome), bbox_inches='tight')
    plt.clf()

    #calculate auprc for each model
    rf_auprc = average_precision_score(rf_y_true, rf_y_pred)
    xgb_auprc = average_precision_score(xgb_y_true, xgb_y_pred)
    keras_auprc = average_precision_score(keras_y_true, keras_y_pred)
    lr_auprc = average_precision_score(lr_y_true, lr_y_pred)

    #calculate auprc 95% ci for each model
    def auprc_ci(y_true, y_pred):
        precision, recall, _ = precision_recall_curve(y_true, y_pred)
        pr_auc = auc(recall, precision)
        mean = pr_auc
        std = sqrt(pr_auc * (1.0 - pr_auc) / len(y_true))
        low  = mean - std
        high = mean + std
        return low, high
    rf_auprc_ci = auprc_ci(rf_y_true, rf_y_pred)
    xgb_auprc_ci = auprc_ci(xgb_y_true, xgb_y_pred)
    keras_auprc_ci = auprc_ci(keras_y_true, keras_y_pred)
    lr_auprc_ci = auprc_ci(lr_y_true, lr_y_pred)

    #calculate precision and recall for each model
    rf_precision, rf_recall, _ = precision_recall_curve(rf_y_true, rf_y_pred)
    xgb_precision, xgb_recall, _ = precision_recall_curve(xgb_y_true, xgb_y_pred)
    keras_precision, keras_recall, _ = precision_recall_curve(keras_y_true, keras_y_pred)
    lr_precision, lr_recall, _ = precision_recall_curve(lr_y_true, lr_y_pred)
    #create labels for precision recall curves
    rf_prc_label = 'RF: ' + str(round(rf_auprc, 3)) + ' (' + str(round(rf_auprc_ci[0], 3)) + ' - ' + str(round(rf_auprc_ci[1], 3)) + ')'
    xgb_prc_label = 'XGB: ' + str(round(xgb_auprc, 3)) + ' (' + str(round(xgb_auprc_ci[0], 3)) + ' - ' + str(round(xgb_auprc_ci[1], 3)) + ')'
    keras_prc_label = 'Keras: ' + str(round(keras_auprc, 3)) + ' (' + str(round(keras_auprc_ci[0], 3)) + ' - ' + str(round(keras_auprc_ci[1], 3)) + ')'
    lr_prc_label = 'LASSO: ' + str(round(lr_auprc, 3)) + ' (' + str(round(lr_auprc_ci[0], 3)) + ' - ' + str(round(lr_auprc_ci[1], 3)) + ')'

    #plot the precision recall curves for each model
    matplotlib.rcParams.update({'font.size': 16})
    #plot the ROC curves for each model
    plt.figure(figsize=(10,10))
    plt.plot(lr_recall, lr_precision, color='cadetblue', label=lr_prc_label)
    plt.plot(rf_recall, rf_precision, color='deepskyblue', label=rf_prc_label)
    plt.plot(xgb_recall, xgb_precision, color='steelblue', label=xgb_prc_label)
    plt.plot(keras_recall, keras_precision, color='dodgerblue', label=keras_prc_label)
    plt.legend(loc="upper right")
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.xlim([0.0, 1.0])
    plt.savefig('results/prc_clin_{}_norect.png'.format(outcome), bbox_inches='tight')

    #save results to csv
    results = pd.DataFrame({'model': ['LR', 'RF', 'XGB', 'NN'],
                        'AUROC': [str(round(lr_confidence[1], 3)), str(round(rf_confidence[1], 3)), str(round(xgb_confidence[1], 3)), str(round(keras_confidence[1], 3))],
                        'AUROC 95% CI': [str(round(lr_confidence[0], 3)) + ' - ' + str(round(lr_confidence[2], 3)), str(round(rf_confidence[0], 3)) + ' - ' + str(round(rf_confidence[2], 3)), str(round(xgb_confidence[0], 3)) + ' - ' + str(round(xgb_confidence[2], 3)), str(round(keras_confidence[0], 3)) + ' - ' + str(round(keras_confidence[2], 3))],
                        'AUPRC': [str(round(lr_auprc, 3)), str(round(rf_auprc, 3)), str(round(xgb_auprc, 3)), str(round(keras_auprc, 3))],
                        'AUPRC 95% CI': [str(round(lr_auprc_ci[0], 3)) + ' - ' + str(round(lr_auprc_ci[1], 3)), str(round(rf_auprc_ci[0], 3)) + ' - ' + str(round(rf_auprc_ci[1], 3)), str(round(xgb_auprc_ci[0], 3)) + ' - ' + str(round(xgb_auprc_ci[1], 3)), str(round(keras_auprc_ci[0], 3)) + ' - ' + str(round(keras_auprc_ci[1], 3))],
                        })
    results.to_csv('results/results_clin_{}_norect.csv'.format(outcome), index=False)
