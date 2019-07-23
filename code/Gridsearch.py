#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:40:07 2019

@author: mw775
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
import rdkit.Chem
import rdkit.Chem.AllChem as Chem
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn import preprocessing
from itertools import product
import pickle

#import data

df_compounds = pd.read_csv('CAMDA-DILI/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv', delimiter=',')


#generate FPs
smis = df_compounds['standardized_smiles'].tolist()
mols = [Chem.MolFromSmiles(smile) for smile in smis]
fps_bit = [Chem.GetMorganFingerprintAsBitVect(mol,2, nBits=2048) for mol in mols]

#store FP bits as input in X_fp_all
X_fp_all = np.empty([920,2048])
for i,j in enumerate(fps_bit):
    for k,l in enumerate(list(j)):
        X_fp_all[i,k] = l
        
# labels (0,1) in Y_all
Y_all = np.empty([920,1])
for i,j in enumerate(df_compounds['vDILIConcern']):
    if j == 'vMost-DILI-Concern' or j == 'vLess-DILI-Concern':
        Y_all[i,0] = 1
        
    if j == 'vNo-DILI-Concern' or j == 'sider_inactive':
        Y_all[i,0] = 0
Y_all = np.squeeze(Y_all)

mc = list(range(174))
lc = list(range(174,434))
nc = list(range(434,661))
sider = list(range(661,920))

#X_fp_mcnc

X_fp_mcnc = np.concatenate((X_fp_all[mc,:],X_fp_all[nc,:]), axis=0)

#Y_mcnc

Y_mcnc = np.concatenate((Y_all[mc],Y_all[nc]))

#X_fp_mclcnc

X_fp_mclcnc = np.concatenate((X_fp_all[mc,:],X_fp_all[lc,:],X_fp_all[nc,:]), axis=0)

#Y_mclcnc

Y_mclcnc = np.concatenate((Y_all[mc],Y_all[lc],Y_all[nc]))


#import descriptors
df_descr = pd.read_csv('CAMDA-DILI/processed_data/Molecular_Descriptors/mol_descriptors_training.csv', delimiter=',', index_col=0)

X_descr_all = df_descr.values
X_descr_all = preprocessing.scale(X_descr_all)

X_descr_mcnc = np.concatenate((X_descr_all[mc,:],X_descr_all[nc,:]), axis=0)

X_descr_mclcnc = np.concatenate((X_descr_all[mc,:],X_descr_all[lc,:],X_descr_all[nc,:]), axis=0)


#import predicted targets
df_targets = pd.read_csv('CAMDA-DILI/processed_data/Target Prediction/training_set_predicted_targets.txt', delimiter='\t')
df_targets = df_targets.set_index('Uniprot')
df_targets = df_targets.T
df_targets = df_targets.iloc[15:,:]
index = []
for i,row in df_targets.iterrows():
    index.append(i.split('\t')[1])
df_targets['new_index'] = index
df_targets = df_targets.set_index('new_index')

dfl = df_compounds.copy()
dfl = dfl.set_index('PubChem_CID')

labels = []
for i,row in df_targets.iterrows():
    labels.append(dfl['vDILIConcern'].loc[int(i)])
    
mc_t = []
lc_t=[]
nc_t =[]
sider_t=[]
for i,label in enumerate(labels): 
    if label == 'vMost-DILI-Concern':
        mc_t.append(i)
    if label == 'vLess-DILI-Concern':
        lc_t.append(i)
    if label == 'vNo-DILI-Concern':
        nc_t.append(i)
    if label == 'sider_inactive':
        sider_t.append(i)
    
X_targets_all = df_targets.values

X_targets_mcnc = np.concatenate((X_targets_all[mc_t],X_targets_all[nc_t]))

X_targets_mclcnc = np.concatenate((X_targets_all[mc_t],X_targets_all[lc_t],X_targets_all[nc_t]))


labels_bin = []
for label in labels:
    if label == 'vMost-DILI-Concern':
        labels_bin.append(1)
    elif label == 'vLess-DILI-Concern':
        labels_bin.append(1)
    else:
        labels_bin.append(0)

Y_targets_all = np.array(labels_bin, dtype='float')
Y_targets_all = np.squeeze(Y_targets_all)

Y_targets_mcnc = np.concatenate((Y_targets_all[mc_t],Y_targets_all[nc_t]))
Y_targets_mclcnc = np.concatenate((Y_targets_all[mc_t],Y_targets_all[lc_t],Y_targets_all[nc_t]))

with open('CAMDA-DILI/processed_data/grid_search/cv_splits/train_splits_pt.pkl', 'rb') as f:
    train_splits_pt_all = pickle.load(f)

with open('CAMDA-DILI/processed_data/grid_search/cv_splits/test_splits_pt.pkl', 'rb') as f:
    test_splits_pt_all = pickle.load(f)

with open('CAMDA-DILI/processed_data/grid_search/cv_splits/train_splits_pt_mcnc.pkl', 'rb') as f:
    train_splits_pt_mcnc = pickle.load(f)

with open('CAMDA-DILI/processed_data/grid_search/cv_splits/test_splits_pt_mcnc.pkl', 'rb') as f:
    test_splits_pt_mcnc = pickle.load(f)

with open('CAMDA-DILI/processed_data/grid_search/cv_splits/train_splits_mclcnc.pkl', 'rb') as f:
    train_splits_mclcnc = pickle.load(f)

with open('CAMDA-DILI/processed_data/grid_search/cv_splits/test_splits_mclcnc.pkl', 'rb') as f:
    test_splits_mclcnc = pickle.load(f)

with open('CAMDA-DILI/processed_data/grid_search/cv_splits/train_splits_pt_mclcnc.pkl', 'rb') as f:
    train_splits_pt_mclcnc = pickle.load(f)

with open('CAMDA-DILI/processed_data/grid_search/cv_splits/test_splits_pt_mclcnc.pkl', 'rb') as f:
    test_splits_pt_mclcnc = pickle.load(f)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state = 1)

# grid_search funtions
def rf_param_selection(X, y, cvsplit=skf):
    est = [100,200,300,400,500,750,1000]
    min_sam = [1,2,3]
    max_dep = [10,15,None]
    param_grid = {'n_estimators': est, 'min_samples_leaf' : min_sam, 'max_depth' : max_dep}
    clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample', bootstrap=True)
    grid_search = GridSearchCV(clf, param_grid, cv=cvsplit, scoring='balanced_accuracy')
    grid_search.fit(X, y)
    grid_search.best_params_
    return(grid_search.best_params_)

def svc_param_selection(X, y, cvsplit=skf):
    Cs = [0.1,0.3, 1,3, 10]
    gammas = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1,'auto']
    param_grid = {'C': Cs, 'gamma' : gammas}
    clf = SVC(class_weight='balanced', random_state=22, kernel='linear')
    grid_search = GridSearchCV(clf, param_grid, cv=cvsplit, scoring='balanced_accuracy')
    grid_search.fit(X, y)
    grid_search.best_params_
    return(grid_search.best_params_)


best_params = {}

methods = [('ECFP','all', 'RF'),('ECFP','mclcnc','RF') ,('ECFP','mcnc','RF'),
           ('ECFP','all','SVM'), ('ECFP','mclcnc','SVM'),
           ('ECFP','mcnc','SVM'),('MD','all','RF'), ('MD','mclcnc','RF'), ('MD','mcnc','RF'),
           ('MD','all','SVM'), ('MD','mclcnc','SVM'), ('MD','mcnc','SVM'),
           ('PT','all','RF'),  ('PT','mclcnc','RF'),
           ('PT','mcnc','RF'),('PT','all','SVM'),('PT','mclcnc','SVM'), ('PT','mcnc','SVM')]

split_count = [1,2,3,4,5]

splits = product(methods,split_count)

bal_acc = []
rec= []
pre = []
auc= []

#'ECFP','all', 'RF'
print(('ECFP','all', 'RF'))    
best_params[('ECFP','all', 'RF')] = rf_param_selection(X_fp_all,Y_all)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('ECFP','all', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('ECFP','all', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('ECFP','all', 'RF')]['max_depth'])

for train, test in skf.split(X_fp_all,Y_all):
    clf.fit(X_fp_all[train],Y_all[train])
    bal_acc.append(balanced_accuracy_score(Y_all[test], clf.predict(X_fp_all[test])))
    rec.append(recall_score(Y_all[test], clf.predict(X_fp_all[test]), pos_label=1))
    pre.append(precision_score(Y_all[test], clf.predict(X_fp_all[test]), pos_label=1))
    auc.append(roc_auc_score(Y_all[test], clf.predict(X_fp_all[test])))

#'ECFP','mclcnc', 'RF'
print(('ECFP','mclcnc', 'RF'))    
best_params[('ECFP','mclcnc', 'RF')] = rf_param_selection(X_fp_mclcnc,Y_mclcnc)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('ECFP','mclcnc', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('ECFP','mclcnc', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('ECFP','mclcnc', 'RF')]['max_depth'])

for train, test in zip(train_splits_mclcnc,test_splits_mclcnc):
    clf.fit(X_fp_mclcnc[train],Y_mclcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_mclcnc[test], clf.predict(X_fp_mclcnc[test])))
    rec.append(recall_score(Y_mclcnc[test], clf.predict(X_fp_mclcnc[test]), pos_label=1))
    pre.append(precision_score(Y_mclcnc[test], clf.predict(X_fp_mclcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_mclcnc[test], clf.predict(X_fp_mclcnc[test])))

#'ECFP','mcnc', 'RF'
print(('ECFP','mcnc', 'RF'))    
best_params[('ECFP','mcnc', 'RF')] = rf_param_selection(X_fp_mcnc,Y_mcnc)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('ECFP','mcnc', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('ECFP','mcnc', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('ECFP','mcnc', 'RF')]['max_depth'])

for train, test in skf.split(X_fp_mcnc,Y_mcnc):
    clf.fit(X_fp_mcnc[train],Y_mcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_mcnc[test], clf.predict(X_fp_mcnc[test])))
    rec.append(recall_score(Y_mcnc[test], clf.predict(X_fp_mcnc[test]), pos_label=1))
    pre.append(precision_score(Y_mcnc[test], clf.predict(X_fp_mcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_mcnc[test], clf.predict(X_fp_mcnc[test])))
    
#'ECFP','all', 'SVM'
print(('ECFP','all', 'SVM'))    
best_params[('ECFP','all', 'SVM')] = svc_param_selection(X_fp_all,Y_all)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('ECFP','all', 'SVM')]['C'],
                                            gamma=best_params[('ECFP','all', 'SVM')]['gamma'])

for train, test in skf.split(X_fp_all,Y_all):
    clf.fit(X_fp_all[train],Y_all[train])
    bal_acc.append(balanced_accuracy_score(Y_all[test], clf.predict(X_fp_all[test])))
    rec.append(recall_score(Y_all[test], clf.predict(X_fp_all[test]), pos_label=1))
    pre.append(precision_score(Y_all[test], clf.predict(X_fp_all[test]), pos_label=1))
    auc.append(roc_auc_score(Y_all[test], clf.predict(X_fp_all[test])))
    
#'ECFP','mclcnc', 'SVM'
print(('ECFP','mclcnc', 'SVM'))    
best_params[('ECFP','mclcnc', 'SVM')] = svc_param_selection(X_fp_mclcnc,Y_mclcnc)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('ECFP','mclcnc', 'SVM')]['C'],
                                            gamma=best_params[('ECFP','mclcnc', 'SVM')]['gamma'])

for train, test in zip(train_splits_mclcnc,test_splits_mclcnc):
    clf.fit(X_fp_mclcnc[train],Y_mclcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_mclcnc[test], clf.predict(X_fp_mclcnc[test])))
    rec.append(recall_score(Y_mclcnc[test], clf.predict(X_fp_mclcnc[test]), pos_label=1))
    pre.append(precision_score(Y_mclcnc[test], clf.predict(X_fp_mclcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_mclcnc[test], clf.predict(X_fp_mclcnc[test])))

#'ECFP','mcnc', 'SVM'
print(('ECFP','mcnc', 'SVM'))    
best_params[('ECFP','mcnc', 'SVM')] = svc_param_selection(X_fp_mcnc,Y_mcnc)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('ECFP','mcnc', 'SVM')]['C'],
                                            gamma=best_params[('ECFP','mcnc', 'SVM')]['gamma'])

for train, test in skf.split(X_fp_mcnc,Y_mcnc):
    clf.fit(X_fp_mcnc[train],Y_mcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_mcnc[test], clf.predict(X_fp_mcnc[test])))
    rec.append(recall_score(Y_mcnc[test], clf.predict(X_fp_mcnc[test]), pos_label=1))
    pre.append(precision_score(Y_mcnc[test], clf.predict(X_fp_mcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_mcnc[test], clf.predict(X_fp_mcnc[test])))
    
    
#'MD','all', 'RF'
print(('MD','all', 'RF'))    
best_params[('MD','all', 'RF')] = rf_param_selection(X_descr_all,Y_all)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('MD','all', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('MD','all', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('MD','all', 'RF')]['max_depth'])

for train, test in skf.split(X_descr_all,Y_all):
    clf.fit(X_descr_all[train],Y_all[train])
    bal_acc.append(balanced_accuracy_score(Y_all[test], clf.predict(X_descr_all[test])))
    rec.append(recall_score(Y_all[test], clf.predict(X_descr_all[test]), pos_label=1))
    pre.append(precision_score(Y_all[test], clf.predict(X_descr_all[test]), pos_label=1))
    auc.append(roc_auc_score(Y_all[test], clf.predict(X_descr_all[test])))
    
#'MD','mclcnc', 'RF'
print(('MD','mclcnc', 'RF'))    
best_params[('MD','mclcnc', 'RF')] = rf_param_selection(X_descr_mclcnc,Y_mclcnc)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('MD','mclcnc', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('MD','mclcnc', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('MD','mclcnc', 'RF')]['max_depth'])

for train, test in zip(train_splits_mclcnc,test_splits_mclcnc):
    clf.fit(X_descr_mclcnc[train],Y_mclcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_mclcnc[test], clf.predict(X_descr_mclcnc[test])))
    rec.append(recall_score(Y_mclcnc[test], clf.predict(X_descr_mclcnc[test]), pos_label=1))
    pre.append(precision_score(Y_mclcnc[test], clf.predict(X_descr_mclcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_mclcnc[test], clf.predict(X_descr_mclcnc[test])))

#'MD','mcnc', 'RF'
print(('MD','mcnc', 'RF'))    
best_params[('MD','mcnc', 'RF')] = rf_param_selection(X_descr_mcnc,Y_mcnc)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('MD','mcnc', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('MD','mcnc', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('MD','mcnc', 'RF')]['max_depth'])

for train, test in skf.split(X_descr_mcnc,Y_mcnc):
    clf.fit(X_descr_mcnc[train],Y_mcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_mcnc[test], clf.predict(X_descr_mcnc[test])))
    rec.append(recall_score(Y_mcnc[test], clf.predict(X_descr_mcnc[test]), pos_label=1))
    pre.append(precision_score(Y_mcnc[test], clf.predict(X_descr_mcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_mcnc[test], clf.predict(X_descr_mcnc[test])))
    

#'MD','all', 'SVM'
print(('MD','all', 'SVM'))    
best_params[('MD','all', 'SVM')] = svc_param_selection(X_descr_all,Y_all)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('MD','all', 'SVM')]['C'],
                                            gamma=best_params[('MD','all', 'SVM')]['gamma'])

for train, test in skf.split(X_descr_all,Y_all):
    clf.fit(X_descr_all[train],Y_all[train])
    bal_acc.append(balanced_accuracy_score(Y_all[test], clf.predict(X_descr_all[test])))
    rec.append(recall_score(Y_all[test], clf.predict(X_descr_all[test]), pos_label=1))
    pre.append(precision_score(Y_all[test], clf.predict(X_descr_all[test]), pos_label=1))
    auc.append(roc_auc_score(Y_all[test], clf.predict(X_descr_all[test])))
    

#'MD','mclcnc', 'SVM'
print(('MD','mclcnc', 'SVM'))    
best_params[('MD','mclcnc', 'SVM')] = svc_param_selection(X_descr_mclcnc,Y_mclcnc)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('MD','mclcnc', 'SVM')]['C'],
                                            gamma=best_params[('MD','mclcnc', 'SVM')]['gamma'])

for train, test in zip(train_splits_mclcnc,test_splits_mclcnc):
    clf.fit(X_descr_mclcnc[train],Y_mclcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_mclcnc[test], clf.predict(X_descr_mclcnc[test])))
    rec.append(recall_score(Y_mclcnc[test], clf.predict(X_descr_mclcnc[test]), pos_label=1))
    pre.append(precision_score(Y_mclcnc[test], clf.predict(X_descr_mclcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_mclcnc[test], clf.predict(X_descr_mclcnc[test])))
    
#'MD','mcnc', 'SVM'
print(('MD','mcnc', 'SVM'))    
best_params[('MD','mcnc', 'SVM')] = svc_param_selection(X_descr_mcnc,Y_mcnc)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('MD','mcnc', 'SVM')]['C'],
                                            gamma=best_params[('MD','mcnc', 'SVM')]['gamma'])

for train, test in skf.split(X_descr_mcnc,Y_mcnc):
    clf.fit(X_descr_mcnc[train],Y_mcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_mcnc[test], clf.predict(X_descr_mcnc[test])))
    rec.append(recall_score(Y_mcnc[test], clf.predict(X_descr_mcnc[test]), pos_label=1))
    pre.append(precision_score(Y_mcnc[test], clf.predict(X_descr_mcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_mcnc[test], clf.predict(X_descr_mcnc[test])))
    

#'PT','all', 'RF'
print(('PT','all', 'RF'))    
best_params[('PT','all', 'RF')] = rf_param_selection(X_targets_all,Y_targets_all)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('PT','all', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('PT','all', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('PT','all', 'RF')]['max_depth'])

for train, test in zip(train_splits_pt_all,test_splits_pt_all):
    clf.fit(X_targets_all[train],Y_targets_all[train])
    bal_acc.append(balanced_accuracy_score(Y_targets_all[test], clf.predict(X_targets_all[test])))
    rec.append(recall_score(Y_targets_all[test], clf.predict(X_targets_all[test]), pos_label=1))
    pre.append(precision_score(Y_targets_all[test], clf.predict(X_targets_all[test]), pos_label=1))
    auc.append(roc_auc_score(Y_targets_all[test], clf.predict(X_targets_all[test])))
    
#'PT','mclcnc', 'RF'
print(('PT','mclcnc', 'RF'))    
best_params[('PT','mclcnc', 'RF')] = rf_param_selection(X_targets_mclcnc,Y_targets_mclcnc)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('PT','mclcnc', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('PT','mclcnc', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('PT','mclcnc', 'RF')]['max_depth'])

for train, test in zip(train_splits_pt_mclcnc,test_splits_pt_mclcnc):
    clf.fit(X_targets_mclcnc[train],Y_targets_mclcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_targets_mclcnc[test], clf.predict(X_targets_mclcnc[test])))
    rec.append(recall_score(Y_targets_mclcnc[test], clf.predict(X_targets_mclcnc[test]), pos_label=1))
    pre.append(precision_score(Y_targets_mclcnc[test], clf.predict(X_targets_mclcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_targets_mclcnc[test], clf.predict(X_targets_mclcnc[test])))

#'PT','mcnc', 'RF'
print(('PT','mcnc', 'RF'))    
best_params[('PT','mcnc', 'RF')] = rf_param_selection(X_targets_mcnc,Y_targets_mcnc)
clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[('PT','mcnc', 'RF')]['n_estimators'],
                              min_samples_leaf=best_params[('PT','mcnc', 'RF')]['min_samples_leaf'],
                              max_depth=best_params[('PT','mcnc', 'RF')]['max_depth'])

for train, test in zip(train_splits_pt_mcnc,test_splits_pt_mcnc):
    clf.fit(X_targets_mcnc[train],Y_targets_mcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_targets_mcnc[test], clf.predict(X_targets_mcnc[test])))
    rec.append(recall_score(Y_targets_mcnc[test], clf.predict(X_targets_mcnc[test]), pos_label=1))
    pre.append(precision_score(Y_targets_mcnc[test], clf.predict(X_targets_mcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_targets_mcnc[test], clf.predict(X_targets_mcnc[test])))
    

#'PT','all', 'SVM'
print(('PT','all', 'SVM'))    
best_params[('PT','all', 'SVM')] = svc_param_selection(X_targets_all,Y_targets_all)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('PT','all', 'SVM')]['C'],
                                            gamma=best_params[('PT','all', 'SVM')]['gamma'])

for train, test in zip(train_splits_pt_all,test_splits_pt_all):
    clf.fit(X_targets_all[train],Y_targets_all[train])
    bal_acc.append(balanced_accuracy_score(Y_targets_all[test], clf.predict(X_targets_all[test])))
    rec.append(recall_score(Y_targets_all[test], clf.predict(X_targets_all[test]), pos_label=1))
    pre.append(precision_score(Y_targets_all[test], clf.predict(X_targets_all[test]), pos_label=1))
    auc.append(roc_auc_score(Y_targets_all[test], clf.predict(X_targets_all[test])))
    
#'PT','mclcnc', 'SVM'
print(('PT','mclcnc', 'SVM'))    
best_params[('PT','mclcnc', 'SVM')] = svc_param_selection(X_targets_mclcnc,Y_targets_mclcnc)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('PT','mclcnc', 'SVM')]['C'],
                                            gamma=best_params[('PT','mclcnc', 'SVM')]['gamma'])

for train, test in zip(train_splits_pt_mclcnc,test_splits_pt_mclcnc):
    clf.fit(X_targets_mclcnc[train],Y_targets_mclcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_targets_mclcnc[test], clf.predict(X_targets_mclcnc[test])))
    rec.append(recall_score(Y_targets_mclcnc[test], clf.predict(X_targets_mclcnc[test]), pos_label=1))
    pre.append(precision_score(Y_targets_mclcnc[test], clf.predict(X_targets_mclcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_targets_mclcnc[test], clf.predict(X_targets_mclcnc[test])))

#'PT','mcnc', 'SVM'
print(('PT','mcnc', 'SVM'))    
best_params[('PT','mcnc', 'SVM')] = svc_param_selection(X_targets_mcnc,Y_targets_mcnc)
clf = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[('PT','mcnc', 'SVM')]['C'],
                                            gamma=best_params[('PT','mcnc', 'SVM')]['gamma'])

for train, test in zip(train_splits_pt_mcnc,test_splits_pt_mcnc):
    clf.fit(X_targets_mcnc[train],Y_targets_mcnc[train])
    bal_acc.append(balanced_accuracy_score(Y_targets_mcnc[test], clf.predict(X_targets_mcnc[test])))
    rec.append(recall_score(Y_targets_mcnc[test], clf.predict(X_targets_mcnc[test]), pos_label=1))
    pre.append(precision_score(Y_targets_mcnc[test], clf.predict(X_targets_mcnc[test]), pos_label=1))
    auc.append(roc_auc_score(Y_targets_mcnc[test], clf.predict(X_targets_mcnc[test])))
    
    
#export cv results
df_cv = pd.DataFrame()
df_cv['split'] = splits
df_cv['balanced_accuracy'] = bal_acc
df_cv['recall'] = rec
df_cv['precision'] = pre
df_cv['ROC_AUC'] = auc

df_cv.to_csv(path_or_buf='CAMDA-DILI/processed_data/grid_search/cv_scores.csv', sep=',')

#export best params
output = open('CAMDA-DILI/processed_data/grid_search/cv_splits/best_params.pkl','wb')

pickle.dump(best_params, output)

output.close
