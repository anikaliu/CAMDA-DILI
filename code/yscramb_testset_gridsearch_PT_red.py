#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 15:07:38 2019

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
df = pd.read_csv('CAMDA-DILI/processed_data/Models/PT/training_set_predicted_targets.txt', delimiter='\t')
df = df.T
df = df.iloc[16:,:]
df.reset_index(inplace=True,drop=True)
df = df.iloc[:659,:]

df_compounds = pd.read_csv('CAMDA-DILI/processed_data/Models/ECFP/standardized_compounds_excl_ambiguous_cluster.csv', delimiter=',',index_col=0)
df_compounds = df_compounds.iloc[:661,:]

#remove compounds not accepted by Pidgin from df compounds and reset index
df_compounds_red = pd.concat([df_compounds.iloc[:135,:],df_compounds.iloc[136:632,:],df_compounds.iloc[633:,:]],axis=0)
df_compounds_red.reset_index(inplace=True,drop=True)

df['vDILIConcern'] = df_compounds_red['vDILIConcern']
df['cluster'] = df_compounds_red['cluster']



#get mappings
#import dictionaries for cluster and compounds
pkl_file = open('CAMDA-DILI/processed_data/Models/GEX/gex_compound_identifier.pkl', 'rb')
compound_ident = pickle.load(pkl_file)
pkl_file.close() 


pkl_file = open('CAMDA-DILI/processed_data/Models/GEX/gex_compound_cluster.pkl', 'rb')
compound_cluster = pickle.load(pkl_file)
pkl_file.close() 

pkl_file = open('CAMDA-DILI/processed_data/Models/GEX/test_list_mcnc.pkl', 'rb')
test_list_mcnc = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('CAMDA-DILI/processed_data/Models/GEX/train_list_mcnc.pkl', 'rb')
train_list_mcnc = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('CAMDA-DILI/processed_data/Models/GEX/test_list_all.pkl', 'rb')
test_list_all = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('CAMDA-DILI/processed_data/Models/GEX/train_list_all.pkl', 'rb')
train_list_all = pickle.load(pkl_file)
pkl_file.close()

#reduce indices in train_mcnc and test _mcnc by 90
for i,split in enumerate(train_list_mcnc):
    for j, comp in enumerate(split):
        train_list_mcnc[i][j]-=90

for i,split in enumerate(test_list_mcnc):
    for j, comp in enumerate(split):
        test_list_mcnc[i][j]-=90


#mapping identifier --> compound
ident_compound = {}
for key in compound_ident:
    ident_compound[compound_ident[key]] = key
#mapping identifier --> cluster 
ident_cluster = {}
for key in ident_compound:
    ident_cluster[key] = compound_cluster[ident_compound[key]]




#mapping compound name -> identifier in chemical data space 
compound_identchem = {}
for i,row in df_compounds.iterrows():
    compound_identchem[row['Compound Name']] = i


  #mapping ident_chem_ident
ident = list(range(178))
ident_chem = []
for i in ident:
    ident_chem.append(compound_identchem[ident_compound[i]]) 

ident_chem_ident = {}
for i,j in zip(ident_chem,ident):
    ident_chem_ident[i] = j

#create mapping all -> all-PT
all_allpt = {}
for i in range(661):
    if i <135:
        all_allpt[i] = i
    if i>135 and i<632:
        all_allpt[i] = i-1
    if i >632:
        all_allpt[i] = i-2
        
ident_chempt = []        
for i in ident_chem:
    ident_chempt.append(all_allpt[i])

#pick subset
df = df.iloc[ident_chempt,:]
df.reset_index(inplace=True,drop=True)


#store FP bits as input in X_fp_all
X_all = df.iloc[:,:-2].to_numpy()

        
# labels (0,1) in Y_all
Y_all = np.empty([178,1])
for i,j in enumerate(df['vDILIConcern']):
    if j == 'vMost-DILI-Concern' or j == 'vLess-DILI-Concern':
        Y_all[i,0] = 1
        
    if j == 'vNo-DILI-Concern':
        Y_all[i,0] = 0
Y_all = np.squeeze(Y_all)

mc = list(range(90,127))
lc = list(range(90))
nc = list(range(127,178))


#X_fp_mcnc

X_mcnc = np.concatenate((X_all[mc,:],X_all[nc,:]), axis=0)

#Y_mcnc

Y_mcnc = np.concatenate((Y_all[mc],Y_all[nc]))

cluster = df['cluster'].values
cluster_mcnc = np.concatenate((cluster[mc],cluster[nc]))


#data for ambis
df_ambis = pd.read_csv('CAMDA-DILI/processed_data/Models/PT/ambiguous_predicted_targets.txt', delimiter='\t')
df_ambis = df_ambis.T
X_ambis = df_ambis.iloc[16:,:].to_numpy()

df_pred = pd.read_csv('CAMDA-DILI/processed_data/Models/myname_predictions_no1_TEMPLATE.txt', delimiter=',')
ambis = df_pred['Compound.Name'].tolist()

# grid_search funtions
def rf_param_selection(X, y, cvsplit):
    est = [100,200,300,400,500,750,1000]
    min_sam = [1,2,3]
    max_dep = [10,15,None]
    param_grid = {'n_estimators': est, 'min_samples_leaf' : min_sam, 'max_depth' : max_dep}
    clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample', bootstrap=True)
    grid_search = GridSearchCV(clf, param_grid, cv=cvsplit, scoring='balanced_accuracy', n_jobs=1)
    grid_search.fit(X, y)
    grid_search.best_params_
    return(grid_search.best_params_)

def svc_param_selection(X, y, cvsplit):
    Cs = [0.05,0.1,0.2,0.3,0.4,0.5, 1]
    param_grid = {'C': Cs}
    clf = SVC(class_weight='balanced', random_state=22, kernel='linear')
    grid_search = GridSearchCV(clf, param_grid, cv=cvsplit, scoring='balanced_accuracy', n_jobs=1)
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

#define random state that will be used for Group-kfold-splits
np.random.seed(55)
state = np.random.get_state()


#0: MCNC, 1: MCLCNC, 2: all
dict_dataset = {0: 'MCNC.',1: 'all.'}
xlist = [X_mcnc,X_all]
ylist = [Y_mcnc,Y_all]
clusterlist = [cluster_mcnc,cluster]
splits_tr = [train_list_mcnc, train_list_all]
splits_te = [test_list_mcnc, test_list_all]



for dataset in range(2):
    X_in = xlist[dataset]
    size = X_in.shape[0]
    Y_in = ylist[dataset]
    cluster_in = clusterlist[dataset]    
    
    #true label and 3 times scrambled

    for i in range(4):
        print('PT','data:',dataset,'scramble',i)
        if i == 0:
            Y = Y_in
        else:
            Y = shuffle(Y_in, random_state = i*20)
            
        
        #skf   
        for j,(train,test) in enumerate(zip(splits_tr[dataset],splits_te[dataset])):
            
                       
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
        
#export predictions in zip file
zf = zipfile.ZipFile('CAMDA-DILI/processed_data/Models/PT_red/Predictions_PT_red.zip', mode='w')
for i,j in zip(predictions,predictions_ident):
    name = j.replace('.','_')
    i.to_csv(path_or_buf='CAMDA-DILI/processed_data/Models/PT_red/Predictions_PT_red/'+name+'.txt',sep=',',index=False)
    zf.write('CAMDA-DILI/processed_data/Models/PT_red/Predictions_PT_red/'+name+'.txt')
    os.remove('CAMDA-DILI/processed_data/Models/PT_red/Predictions_PT_red/'+name+'.txt')
zf.close()
    
df_ts = pd.DataFrame()
df_ts['splits'] = ts_splits
df_ts['accuracy'] = acc_ts
df_ts['balanced_accuracy'] = bal_acc_ts
df_ts['recall'] = rec_ts
df_ts['precision'] = pre_ts
df_ts['AU_Prec_Rec_Curve'] = aupr_ts
df_ts['ROC_AUC'] = rocauc_ts
df_ts['MCC'] = mcc_ts


df_ts.to_csv('CAMDA-DILI/processed_data/Models/PT_red/ts_scores_PT_red_yscr.csv',sep=',',index=False)

df_cv = pd.DataFrame()
df_cv['splits'] = cv_splits
df_cv['accuracy'] = acc
df_cv['balanced_accuracy'] = bal_acc
df_cv['recall'] = rec
df_cv['precision'] = pre
df_cv['AU_Prec_Rec_Curve'] = aupr
df_cv['ROC_AUC'] = rocauc
df_cv['MCC'] = mcc

df_cv.to_csv('CAMDA-DILI/processed_data/Models/PT_red/cv_scores_PT_red_yscr.csv',sep=',',index=False)

#export best params
output = open('CAMDA-DILI/processed_data/Models/PT_red/best_params_PT_red_yscr.pkl','wb')

pickle.dump(best_params, output)

output.close()
