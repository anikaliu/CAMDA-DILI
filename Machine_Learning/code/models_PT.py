#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 12:23:32 2019

@author: mw775
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import matthews_corrcoef
from sklearn.utils import shuffle
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import GroupKFold
import pickle
import os
import zipfile


#import data
df_compounds = pd.read_csv('CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_excl_ambiguous_cluster.csv', delimiter=',')

df = pd.read_csv('CAMDA-DILI/Data_Processing/Target_Prediction/data/training_set_predicted_targets.txt', delimiter='\t')
df = df.T
X_all = df.iloc[16:,:].to_numpy()


# labels (0,1) in Y_all
Y_all = np.empty([923,1])
for i,j in enumerate(df_compounds['vDILIConcern']):
    if j == 'vMost-DILI-Concern' or j == 'vLess-DILI-Concern':
        Y_all[i,0] = 1
        
    if j == 'vNo-DILI-Concern' or j == 'sider_inactive':
        Y_all[i,0] = 0
Y_all = np.squeeze(Y_all)


#remove compounds not accepted by Pidgin from Y_all
Y_all = np.delete(Y_all,[135,632,665,670,873,922])


#ranges specific for PT
mc = list(range(173))
lc = list(range(173,433))
nc = list(range(433,659))
sider = list(range(659,917))

#X_pt_mcnc

X_mcnc = np.concatenate((X_all[mc,:],X_all[nc,:]), axis=0)

#Y_mcnc

Y_mcnc = np.concatenate((Y_all[mc],Y_all[nc]))

cluster = np.squeeze(df_compounds['cluster'].values)

#remove compounds not accepted by Pidgin from cluster
cluster = np.delete(cluster,[135,632,665,670,873,922])

cluster_mcnc = np.concatenate((cluster[mc],cluster[nc]))

#X&Y_mclcnc
X_mclcnc = np.concatenate((X_all[mc,:],X_all[lc,:],X_all[nc,:]), axis=0)
Y_mclcnc = np.concatenate((Y_all[mc],Y_all[lc],Y_all[nc]))
cluster_mclcnc = np.concatenate((cluster[mc],cluster[lc],cluster[nc]))



df_pred = pd.read_csv('CAMDA-DILI/Data_Processing/Challenge_Data/myname_predictions_no1_TEMPLATE.txt', delimiter=',')

#data for ambis
df_ambis = pd.read_csv('CAMDA-DILI/Data_Processing/Target_Prediction/data/ambiguous_predicted_targets.txt', delimiter='\t')
df_ambis = df_ambis.T
X_ambis = df_ambis.iloc[16:,:].to_numpy()



# grid_search funtions
def rf_param_selection(X, y, cvsplit):
    est = [100,200,300,400,500,750,1000]
    min_sam = [1,2,3]
    max_dep = [10,15,None]
    param_grid = {'n_estimators': est, 'min_samples_leaf' : min_sam, 'max_depth' : max_dep}
    clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample', bootstrap=True)
    grid_search = GridSearchCV(clf, param_grid, cv=cvsplit, scoring='balanced_accuracy', n_jobs=15)
    grid_search.fit(X, y)
    grid_search.best_params_
    return(grid_search.best_params_)

def svc_param_selection(X, y, cvsplit):
    Cs = [0.05,0.1,0.2,0.3,0.4,0.5, 1]
    param_grid = {'C': Cs}
    clf = SVC(class_weight='balanced', random_state=22, kernel='linear')
    grid_search = GridSearchCV(clf, param_grid, cv=cvsplit, scoring='balanced_accuracy', n_jobs=15)
    grid_search.fit(X, y)
    grid_search.best_params_
    return(grid_search.best_params_)
    

best_params = {}

acc = []
bal_acc = []
rec = []
pre = []
rocauc = []
aupr = []
mcc = []


acc_ts = []
bal_acc_ts = []
rec_ts = []
pre_ts = []
rocauc_ts = []
aupr_ts = []
mcc_ts = []

ts_splits = []

cv_splits = []

predictions = []
predictions_ident = []

features_rf = []
features_ident = []

features_svm = []


#define random state that will be used for Group-kfold-splits
np.random.seed(55)
state = np.random.get_state()


#import skf splits
pkl_file = open('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/ECFP/splits/all_ttsplits_tr.pkl', 'rb')
all_ttsplits_tr = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/ECFP/splits/all_ttsplits_te.pkl', 'rb')
all_ttsplits_te = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/ECFP/splits/mclcnc_ttsplits_tr.pkl', 'rb')
mclcnc_ttsplits_tr = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/ECFP/splits/mclcnc_ttsplits_te.pkl', 'rb')
mclcnc_ttsplits_te = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/ECFP/splits/mcnc_ttsplits_tr.pkl', 'rb')
mcnc_ttsplits_tr = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/ECFP/splits/mcnc_ttsplits_te.pkl', 'rb')
mcnc_ttsplits_te = pickle.load(pkl_file)
pkl_file.close()



#create mapping indices all -> all-PT, since 6 compounds not in PT data set (dropped by Pidgin)

all_allpt = {}
for i in range(923):
    if i <135:
        #no change in this range
        all_allpt[i] = i
    if i>135 and i<632:
        all_allpt[i] = i-1
    if i >632 and i < 665:
        all_allpt[i] = i-2
    if i > 665 and i < 670:
        all_allpt[i] = i-3
    if i >670 and i<873:
        all_allpt[i] = i-4
    if i>873 and i<922:
        all_allpt[i] = i-5
        
#mapping mclcnc-> mclcnc-PT
mclcnc_mclcncpt = all_allpt
#create mapping mcnc -> mcnc-PT
mcnc_mcncpt = {}
for i in range(401):
    if i<135:
        mcnc_mcncpt[i] = i
    if i>135 and i<(632-260):
        mcnc_mcncpt[i] = i-1
    if i>(632-260):
        mcnc_mcncpt[i] = i-2

  

#map splits to PT indices
def map_splits(splits_in,mapping_dict):
    new_splits = []
    for i,sp in enumerate(splits_in):
        new_splits.append([])
        for comp in sp:
            if comp in mapping_dict:
                new_splits[i].append(mapping_dict[comp])
            if comp not in mapping_dict:
                continue
    for i,sp in enumerate(new_splits):
        new_splits[i] = np.array(sp)
    return(new_splits)

all_ttsplits_tr = map_splits(all_ttsplits_tr,all_allpt)
all_ttsplits_te = map_splits(all_ttsplits_te,all_allpt)

mclcnc_ttsplits_tr = map_splits(mclcnc_ttsplits_tr,mclcnc_mclcncpt)
mclcnc_ttsplits_te = map_splits(mclcnc_ttsplits_te,mclcnc_mclcncpt)

mcnc_ttsplits_tr = map_splits(mcnc_ttsplits_tr,mcnc_mcncpt)
mcnc_ttsplits_te = map_splits(mcnc_ttsplits_te,mcnc_mcncpt)        



#0: MCNC, 1: MCLCNC, 2: all
dict_dataset = {0: 'MCNC.',1: 'MCLCNC.',2: 'all.'}
xlist = [X_mcnc,X_mclcnc,X_all]
ylist = [Y_mcnc,Y_mclcnc,Y_all]
clusterlist = [cluster_mcnc,cluster_mclcnc,cluster]


for dataset in range(3):
    X_in = xlist[dataset]
    size = X_in.shape[0]
    Y_in = ylist[dataset]
    cluster_in = clusterlist[dataset]    
    splits_tr = [mcnc_ttsplits_tr, mclcnc_ttsplits_tr, all_ttsplits_tr]
    splits_te = [mcnc_ttsplits_te, mclcnc_ttsplits_te, all_ttsplits_te]

    #true label and 3 times scrambled

    for i in range(4):
        print('PT','data:',dataset,'scramble',i)
        if i == 0:
            Y = Y_in
        else:
            Y = shuffle(Y_in, random_state = i*20)
            
        
        #skf   
        for j,(train,test) in enumerate(zip(splits_tr[dataset],splits_te[dataset])):
            
        
            #size of test set: 10%
            test_idx = test
            train_idx = train
    
            #set random state for gkf
            np.random.set_state(state)
            gkf = GroupKFold(n_splits=5).split(X_in[train_idx],Y[train_idx],cluster_in[train_idx])
    
      
            best_params[dict_dataset[dataset]+'RF.' + str(i+1)+'.'+str(j+1)] = rf_param_selection(X_in[train_idx],Y[train_idx],gkf)
    
            np.random.set_state(state)
            gkf = GroupKFold(n_splits=5).split(X_in[train_idx],Y[train_idx],cluster_in[train_idx])
    
    
            best_params[dict_dataset[dataset]+'SVM.' + str(i+1)+'.'+str(j+1)] = svc_param_selection(X_in[train_idx],Y[train_idx],gkf)
    
  
    
            clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[dict_dataset[dataset]+'RF.'+str(i+1)+'.'+str(j+1)]['n_estimators'],
                              min_samples_leaf=best_params[dict_dataset[dataset]+'RF.'+str(i+1)+'.'+str(j+1)]['min_samples_leaf'],
                              max_depth=best_params[dict_dataset[dataset]+'RF.'+str(i+1)+'.'+str(j+1)]['max_depth'])
            clf2 = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[dict_dataset[dataset]+'SVM.'+str(i+1)+'.'+str(j+1)]['C'],
               probability=True)
        
            np.random.set_state(state)
            gkf = GroupKFold(n_splits=5).split(X_in[train_idx],Y[train_idx],cluster_in[train_idx])
            for k, (train,test) in enumerate(gkf):
                clf.fit(X_in[train_idx][train],Y[train_idx][train])
                acc.append(accuracy_score(Y[train_idx][test], clf.predict(X_in[train_idx][test])))
                bal_acc.append(balanced_accuracy_score(Y[train_idx][test], clf.predict(X_in[train_idx][test])))
                rec.append(recall_score(Y[train_idx][test], clf.predict(X_in[train_idx][test]), pos_label=1))
                pre.append(precision_score(Y[train_idx][test], clf.predict(X_in[train_idx][test]), pos_label=1))
                rocauc.append(roc_auc_score(Y[train_idx][test], clf.predict_proba(X_in[train_idx][test])[:,1]))
                precision,recall,_ = precision_recall_curve(Y[train_idx][test], clf.predict_proba(X_in[train_idx][test])[:,1])
                aupr.append(auc(recall,precision))
                mcc.append(matthews_corrcoef(Y[train_idx][test], clf.predict(X_in[train_idx][test])))
                cv_splits.append(dict_dataset[dataset]+'RF.'+str(i+1)+'.'+str(j+1)+'.'+str(k+1))
            
            
                clf2.fit(X_in[train_idx][train],Y[train_idx][train])
                acc.append(accuracy_score(Y[train_idx][test], clf2.predict(X_in[train_idx][test])))
                bal_acc.append(balanced_accuracy_score(Y[train_idx][test], clf2.predict(X_in[train_idx][test])))
                rec.append(recall_score(Y[train_idx][test], clf2.predict(X_in[train_idx][test]), pos_label=1))
                pre.append(precision_score(Y[train_idx][test], clf2.predict(X_in[train_idx][test]), pos_label=1))
                rocauc.append(roc_auc_score(Y[train_idx][test], clf2.predict_proba(X_in[train_idx][test])[:,1]))
                precision,recall,_ = precision_recall_curve(Y[train_idx][test], clf2.predict_proba(X_in[train_idx][test])[:,1])
                aupr.append(auc(recall,precision))
                mcc.append(matthews_corrcoef(Y[train_idx][test], clf2.predict(X_in[train_idx][test])))
                cv_splits.append(dict_dataset[dataset]+'SVM.'+str(i+1)+'.'+str(j+1)+'.'+str(k+1))
                
                
                
            #evaluate on test set
            clf.fit(X_in[train_idx],Y[train_idx])
            acc_ts.append(accuracy_score(Y[test_idx], clf.predict(X_in[test_idx])))
            bal_acc_ts.append(balanced_accuracy_score(Y[test_idx], clf.predict(X_in[test_idx])))
            rec_ts.append(recall_score(Y[test_idx], clf.predict(X_in[test_idx]), pos_label=1))
            pre_ts.append(precision_score(Y[test_idx], clf.predict(X_in[test_idx]), pos_label=1))
            rocauc_ts.append(roc_auc_score(Y[test_idx], clf.predict_proba(X_in[test_idx])[:,1]))
            precision,recall,_ = precision_recall_curve(Y[test_idx], clf.predict_proba(X_in[test_idx])[:,1])
            aupr_ts.append(auc(recall,precision))
            mcc_ts.append(matthews_corrcoef(Y[test_idx], clf.predict(X_in[test_idx])))
            ts_splits.append(dict_dataset[dataset]+'RF.'+str(i+1)+'.'+str(j+1))
            
                        
            clf2.fit(X_in[train_idx],Y[train_idx])
            acc_ts.append(accuracy_score(Y[test_idx], clf2.predict(X_in[test_idx])))
            bal_acc_ts.append(balanced_accuracy_score(Y[test_idx], clf2.predict(X_in[test_idx])))
            rec_ts.append(recall_score(Y[test_idx], clf2.predict(X_in[test_idx]), pos_label=1))
            pre_ts.append(precision_score(Y[test_idx], clf2.predict(X_in[test_idx]), pos_label=1))
            rocauc_ts.append(roc_auc_score(Y[test_idx], clf2.predict_proba(X_in[test_idx])[:,1]))
            precision,recall,_ = precision_recall_curve(Y[test_idx], clf2.predict_proba(X_in[test_idx])[:,1])
            aupr_ts.append(auc(recall,precision))
            mcc_ts.append(matthews_corrcoef(Y[test_idx], clf2.predict(X_in[test_idx])))
            ts_splits.append(dict_dataset[dataset]+'SVM.'+str(i+1)+'.'+str(j+1))

        
            #make predictions for each set with the corresponding parameters
            if i ==0:
                clf.fit(X_in,Y)
                pred = clf.predict(X_ambis)
                df_curr_pred = df_pred.copy()
                for idx,p in enumerate(pred):
                    if p==1:
                        df_curr_pred.iloc[idx,1] = 'DILI'
                    else:
                        df_curr_pred.iloc[idx,1] = 'NoDILI'
                predictions.append(df_curr_pred)
                predictions_ident.append(dict_dataset[dataset]+'RF.'+str(i+1)+'.'+str(j+1))
                
                #get feature importances
                features_rf.append(list(clf.feature_importances_))
                features_ident.append(dict_dataset[dataset]+str(i+1)+'.'+str(j+1))
                
                clf2.fit(X_in,Y)
                pred = clf2.predict(X_ambis)
                df_curr_pred = df_pred.copy()
                for idx,p in enumerate(pred):
                    if p==1:
                        df_curr_pred.iloc[idx,1] = 'DILI'
                    else:
                        df_curr_pred.iloc[idx,1] = 'NoDILI'
                predictions.append(df_curr_pred)
                predictions_ident.append(dict_dataset[dataset]+'SVM.'+str(i+1)+'.'+str(j+1))
                
                features_svm.append(list(clf2.coef_[0]))
                
          
#export predictions in zip file
zf = zipfile.ZipFile('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/Predictions_PT.zip', mode='w')
for i,j in zip(predictions,predictions_ident):
    name = j.replace('.','_')
    i.to_csv(path_or_buf='CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/Predictions_PT/'+name+'.txt',sep=',',index=False)
    zf.write('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/Predictions_PT/'+name+'.txt')
    os.remove('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/Predictions_PT/'+name+'.txt')
zf.close()
    

#export feature_importances
df_feat_rf = pd.DataFrame()
for feat,ident in zip(features_rf,features_ident):
    df_feat_rf[ident] = feat
df_feat_rf.to_csv('CAMDA-DILI/processed_data/Models/PT/CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/rf_feature_importance.csv',sep=',',index=False)

df_feat_svm = pd.DataFrame()
for feat,ident in zip(features_svm,features_ident):
    df_feat_svm[ident] = feat
df_feat_svm.to_csv('CAMDA-DILI/processed_data/Models/PT/CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/svm_feature_coefficients.csv',sep=',',index=False)
   

df_ts = pd.DataFrame()
df_ts['splits'] = ts_splits
df_ts['accuracy'] = acc_ts
df_ts['balanced_accuracy'] = bal_acc_ts
df_ts['recall'] = rec_ts
df_ts['precision'] = pre_ts
df_ts['AU_Prec_Rec_Curve'] = aupr_ts
df_ts['ROC_AUC'] = rocauc_ts
df_ts['MCC'] = mcc_ts


df_ts.to_csv('CAMDA-DILI/processed_data/Models/PT/CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/test_scores_PT.csv',sep=',',index=False)

df_cv = pd.DataFrame()
df_cv['splits'] = cv_splits
df_cv['accuracy'] = acc
df_cv['balanced_accuracy'] = bal_acc
df_cv['recall'] = rec
df_cv['precision'] = pre
df_cv['AU_Prec_Rec_Curve'] = aupr
df_cv['ROC_AUC'] = rocauc
df_cv['MCC'] = mcc

df_cv.to_csv('CAMDA-DILI/processed_data/Models/PT/CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/cv_scores_PT.csv',sep=',',index=False)

#export best params
output = open('CAMDA-DILI/processed_data/Models/PT/CAMDA-DILI/Machine_Learning/data/Model_Results_Parameters/PT/best_parameters_PT.pkl','wb')

pickle.dump(best_params, output)

output.close()

