import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.utils import shuffle
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import matthews_corrcoef
import pickle
import zipfile
import os

#define random state that will be used for Group-kfold-splits
np.random.seed(55)
state = np.random.get_state()


#import dictionariy for mapping from compounds name to unique GEX compound identifier(0-177)
pkl_file = open('CAMDA-DILI/Data_Processing/Gene_Expression/data/gex_compound_identifier.pkl', 'rb')

compound_ident = pickle.load(pkl_file)

pkl_file.close() 


#indices for compound identifier lc: 0-89, mc: 90-126, nc: 127-177
#lc = np.arange(0,90)
mc = np.arange(90,127)
nc = np.arange(127,178)

#create corresponding labels
Y_all = [1 if i<=126 else 0 for i in range(178)]
Y_all = np.array(Y_all)
Y_mcnc = Y_all[90:]

#generate train-test-splits on the level of identifier(0-177), seperately for all data and mcnc
#skf for training - test splitting
train_list_all = []
test_list_all = []
train_list_mcnc = []
test_list_mcnc = []



skf = StratifiedKFold(n_splits=10,shuffle=True)


#mcnc    
np.random.seed(66)
for train,test in skf.split(np.zeros(88),Y_mcnc):
    train_list_mcnc.append(train)
    test_list_mcnc.append(test)
    
#add 90 to each value of mcnc to get original identifier
for a,li in enumerate(train_list_mcnc):
    for j,k in enumerate(li):
        train_list_mcnc[a][j] = k+90        

for a,li in enumerate(test_list_mcnc):
    for j,k in enumerate(li):
        test_list_mcnc[a][j] = k+90
        

# grid_search funtions
def rf_param_selection(X, y, cvsplit):
    est = [100,200,300,400,500,750,1000]
    min_sam = [1,2,3]
    max_dep = [10,15,None]
    param_grid = {'n_estimators': est, 'min_samples_leaf' : min_sam, 'max_depth' : max_dep}
    clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample', bootstrap=True)
    grid_search = GridSearchCV(clf, param_grid, cv=cvsplit, scoring='balanced_accuracy', n_jobs=31)
    grid_search.fit(X, y)
    grid_search.best_params_
    return(grid_search.best_params_)

def svc_param_selection(X, y, cvsplit):
    Cs = [0.05,0.1,0.2,0.3,0.4,0.5, 1]
    param_grid = {'C': Cs}
    clf = SVC(class_weight='balanced', random_state=22, kernel='linear')
    grid_search = GridSearchCV(clf, param_grid, cv=cvsplit, scoring='balanced_accuracy', n_jobs=31)
    grid_search.fit(X, y)
    grid_search.best_params_
    return(grid_search.best_params_)
    
#dictinary to save best parameters
best_params = {}

#lists to save performance in cv and on testsets (_ts)
acc = []
bal_acc = []
rec = []
pre = []
rocauc = []
aupr = []
mcc = []
cv_exceptions = []
cv_exception_ytrue = []

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

#import data
basedir = 'data_files'

files = ["('A375', 6, 10.0).csv", "('A549', 24, 10.0).csv", "All_Cell_Lines.csv", "('ASC', 24, 10.0).csv", "('HA1E', 6, 10.0).csv", "('HCC515', 6, 10.0).csv", "('HEPG2', 6, 10.0).csv",
         "('HT29', 6, 10.0).csv", "('MCF7', 6, 10.0).csv", "('PC3', 6, 10.0).csv", "('PC3', 24, 10.0).csv", "('PHH', 24, 10.0).csv", "('SKB', 24, 10.0).csv", "('VCAP', 6, 10.0).csv"]

 
basedir2 = 'data_files_external'


files2 = ["External_('A375', 6, 10.0).csv", "External_('A549', 24, 10.0).csv", "External_All_Cell_Lines.csv", "External_('ASC', 24, 10.0).csv", "External_('HA1E', 6, 10.0).csv",
         "External_('HCC515', 6, 10.0).csv", "External_('HEPG2', 6, 10.0).csv", "External_('HT29', 6, 10.0).csv", "External_('MCF7', 6, 10.0).csv", "External_('MCF7', 24, 10.0).csv", "External_('PC3', 6, 10.0).csv",
         "External_('PC3', 24, 10.0).csv", "External_('PHH', 24, 10.0).csv", "External_('SKB', 24, 10.0).csv", "External_('VCAP', 6, 10.0).csv"]

#template for ambiguous predictions
df_pred_temp = pd.read_csv('myname_predictions_no1_TEMPLATE.txt', delimiter=',')

for file,file_ext in zip(files,files2):
    df = pd.read_csv(basedir+'/'+file,delimiter=',',index_col=[0])
    #df = pd.concat([df[df['vDILIConcern'] == 'vMost-DILI-Concern'], df[df['vDILIConcern'] == 'vNo-DILI-Concern']], ignore_index=True)
    df.drop(['index'],axis=1,inplace=True)
    df = df.sort_values(by=['vDILIConcern','pert_name'])
    #assign compound identifier for each file in a general manner
    #identifier: sequence of identifiers (0-177) in current file with its replicates
    identifier = []
    
    for _,row in df.iterrows():
        identifier.append(compound_ident[row['pert_name']])
    
    labels = []
    for _,row in df.iterrows():
        if row['vDILIConcern'] == 'vLess-DILI-Concern' or row['vDILIConcern'] == 'vMost-DILI-Concern':
            labels.append(1)
        else:
            labels.append(0)
    
    df['identifier'] = identifier
    df.reset_index(drop=True,inplace=True)
    
    print(df.shape)
    #get first appearance of Most Concern
    for idx,row in df.iterrows():
        if row['vDILIConcern'] == 'vMost-DILI-Concern':
            start_mc = idx
            break
    
    #prepare ambiguous
    df_ambi = pd.read_csv(basedir2+'/'+file_ext,delimiter=',',index_col=[0])
    df_ambi.reset_index(drop=True,inplace=True)
    df_ambi = df_ambi.sort_values(by=['pert_name'])
    
    ambi_list = df_ambi['pert_name'].tolist()
    X_ambi = df_ambi.iloc[:,:-2].to_numpy()
    
    
    #iterate through data combinations. #0: all, 1: mcnc
    dict_dataset = {0: 'MCNC.'}
    for dataset in range(1):
        
                   
        if dataset == 0:
            identifier = np.array(identifier[start_mc:])
            X = df.iloc[:,:-4].to_numpy()[start_mc:,:]
            Yin = np.array(labels)[start_mc:]
            #translate identifier in folds into np.array indices
            train_list_curr = []
            test_list_curr = []
            for i,(train,test) in enumerate(zip(train_list_mcnc,test_list_mcnc)):
                train_list_curr.append([])
                test_list_curr.append([])
                for comp_tr in train:
                    indices = [j for j,ident in enumerate(identifier) if ident==comp_tr]
                    for idx in indices:
                        train_list_curr[i].append(idx)
                for comp_te in test:
                    indices = [j for j,ident in enumerate(identifier) if ident==comp_te]
                    for idx in indices:
                        test_list_curr[i].append(idx)
                        
          
        
        
        #yscrambling scr=0: correct labels
        for scr in range(4):
            #print("Scr: ", scr) ## 0505
            #print("Yin: ", Yin) ## 0505
            print(file,dataset,scr)
            if scr == 0:
                Y = Yin
            else:
                Y = shuffle(Yin, random_state = scr*14)
            
            #iterate through different train test splits
            for split,(train, test) in enumerate(zip(train_list_curr,test_list_curr)):
                
                test_idx = test
                train_idx = train
                
                               
                #set random state for gkf
                np.random.set_state(state)
                #print(X[train_idx])
                gkf = GroupKFold(n_splits=5).split(X[train_idx],Y[train_idx],identifier[train_idx])
                
                #get RF parameters
                best_k = best_params[file+'.'+dict_dataset[dataset]+'RF.' + str(scr+1)+'.'+str(split+1)] = rf_param_selection(X[train_idx],Y[train_idx],
                            gkf)
                
                np.random.set_state(state)
                gkf = GroupKFold(n_splits=5).split(X[train_idx],Y[train_idx],identifier[train_idx])
                
                #get SVM parameters
                best_params[file+'.'+dict_dataset[dataset]+'SVM.' + str(scr+1)+'.'+str(split+1)] = svc_param_selection(X[train_idx],Y[train_idx],
                            gkf)
                        
                #clf: RF
                clf = RandomForestClassifier(random_state = 2,class_weight = 'balanced_subsample',
                             bootstrap=True, n_estimators=best_params[file+'.'+dict_dataset[dataset]+'RF.' + str(scr+1)+'.'+str(split+1)]['n_estimators'],
                              min_samples_leaf=best_params[file+'.'+dict_dataset[dataset]+'RF.' + str(scr+1)+'.'+str(split+1)]['min_samples_leaf'],
                              max_depth=best_params[file+'.'+dict_dataset[dataset]+'RF.' + str(scr+1)+'.'+str(split+1)]['max_depth'])
                #clf2: SVM
                clf2 = SVC(class_weight='balanced', random_state=22, kernel='linear', C=best_params[file+'.'+dict_dataset[dataset]+'SVM.' + str(scr+1)+'.'+str(split+1)]['C'],
               probability=True)
                
                
        
                np.random.set_state(state)
                gkf = GroupKFold(n_splits=5).split(X[train_idx],Y[train_idx],identifier[train_idx])
                #get RF/SVM cv performances for best parameter set set 
                for k, (traincv,testcv) in enumerate(gkf):
                    clf.fit(X[train_idx][traincv],Y[train_idx][traincv])
                    try:
                        acc.append(accuracy_score(Y[train_idx][testcv], clf.predict(X[train_idx][testcv])))
                        bal_acc.append(balanced_accuracy_score(Y[train_idx][testcv], clf.predict(X[train_idx][testcv])))
                        rec.append(recall_score(Y[train_idx][testcv], clf.predict(X[train_idx][testcv]), pos_label=1))
                        pre.append(precision_score(Y[train_idx][testcv], clf.predict(X[train_idx][testcv]), pos_label=1))
                        rocauc.append(roc_auc_score(Y[train_idx][testcv], clf.predict_proba(X[train_idx][testcv])[:,1]))
                        precision,recall,_ = precision_recall_curve(Y[train_idx][testcv], clf.predict_proba(X[train_idx][testcv])[:,1])
                        aupr.append(auc(recall,precision))
                        mcc.append(matthews_corrcoef(Y[train_idx][testcv], clf.predict(X[train_idx][testcv])))
                        cv_splits.append(file+'.'+dict_dataset[dataset]+'RF.'+str(scr+1)+'.'+str(split+1)+'.'+str(k+1))
                    #raise exception if only 1 class present, ROCAUC would raise an error then
                    except:
                        acc.append('exception')
                        bal_acc.append('exception')
                        rec.append('exception')
                        pre.append('exception')
                        rocauc.append('exception')
                        aupr.append('exception')
                        mcc.append('exception')
                        cv_splits.append(file+'.'+dict_dataset[dataset]+'RF.'+str(scr+1)+'.'+str(split+1)+'.'+str(k+1))
                        cv_exceptions.append(file+'.'+dict_dataset[dataset]+'RF.'+str(scr+1)+'.'+str(split+1)+'.'+str(k+1))
                        cv_exception_ytrue.append(Y[train_idx][testcv])
                        
                
                    clf2.fit(X[train_idx][traincv],Y[train_idx][traincv])
                    try:
                        acc.append(accuracy_score(Y[train_idx][testcv], clf2.predict(X[train_idx][testcv])))
                        bal_acc.append(balanced_accuracy_score(Y[train_idx][testcv], clf2.predict(X[train_idx][testcv])))
                        rec.append(recall_score(Y[train_idx][testcv], clf2.predict(X[train_idx][testcv]), pos_label=1))
                        pre.append(precision_score(Y[train_idx][testcv], clf2.predict(X[train_idx][testcv]), pos_label=1))
                        rocauc.append(roc_auc_score(Y[train_idx][testcv], clf2.predict_proba(X[train_idx][testcv])[:,1]))
                        precision,recall,_ = precision_recall_curve(Y[train_idx][testcv], clf2.predict_proba(X[train_idx][testcv])[:,1])
                        aupr.append(auc(recall,precision))
                        mcc.append(matthews_corrcoef(Y[train_idx][testcv], clf2.predict(X[train_idx][testcv])))
                        cv_splits.append(file+'.'+dict_dataset[dataset]+'SVM.'+str(scr+1)+'.'+str(split+1)+'.'+str(k+1))
                    except:
                        acc.append('exception')
                        bal_acc.append('exception')
                        rec.append('exception')
                        pre.append('exception')
                        rocauc.append('exception')
                        aupr.append('exception')
                        mcc.append('exception')
                        cv_splits.append(file+'.'+dict_dataset[dataset]+'SVM.'+str(scr+1)+'.'+str(split+1)+'.'+str(k+1))
                        cv_exceptions.append(file+'.'+dict_dataset[dataset]+'SVM.'+str(scr+1)+'.'+str(split+1)+'.'+str(k+1))
                        cv_exception_ytrue.append(Y[train_idx][testcv])        
                
                #get test set performances
                clf.fit(X[train_idx],Y[train_idx])
                acc_ts.append(accuracy_score(Y[test_idx], clf.predict(X[test_idx])))
                bal_acc_ts.append(balanced_accuracy_score(Y[test_idx], clf.predict(X[test_idx])))
                rec_ts.append(recall_score(Y[test_idx], clf.predict(X[test_idx]), pos_label=1))
                pre_ts.append(precision_score(Y[test_idx], clf.predict(X[test_idx]), pos_label=1))
                #rocauc_ts.append(roc_auc_score(Y[test
                
                try:
                    rocauc_ts.append(roc_auc_score(Y[test_idx], clf.predict_proba(X[test_idx])[:,1]))
                except:
                    rocauc_ts.append('exception')

                precision,recall,_ = precision_recall_curve(Y[test_idx], clf.predict_proba(X[test_idx])[:,1])
                aupr_ts.append(auc(recall,precision))
                mcc_ts.append(matthews_corrcoef(Y[test_idx], clf.predict(X[test_idx])))
                ts_splits.append(file+'.'+dict_dataset[dataset]+'RF.'+str(scr+1)+'.'+str(split+1))
                
                
                clf2.fit(X[train_idx],Y[train_idx])
                acc_ts.append(accuracy_score(Y[test_idx], clf2.predict(X[test_idx])))
                bal_acc_ts.append(balanced_accuracy_score(Y[test_idx], clf2.predict(X[test_idx])))
                rec_ts.append(recall_score(Y[test_idx], clf2.predict(X[test_idx]), pos_label=1))
                pre_ts.append(precision_score(Y[test_idx], clf2.predict(X[test_idx]), pos_label=1))
                try:
                    rocauc_ts.append(roc_auc_score(Y[test_idx], clf2.predict_proba(X[test_idx])[:,1]))
                except:
                    rocauc_ts.append('Exception')
                #rocauc_ts.append(roc_auc_score(Y[test_idx], clf2.predict_proba(X[test_idx])[:,1]))
                precision,recall,_ = precision_recall_curve(Y[test_idx], clf2.predict_proba(X[test_idx])[:,1])
                aupr_ts.append(auc(recall,precision))
                mcc_ts.append(matthews_corrcoef(Y[test_idx], clf2.predict(X[test_idx])))
                ts_splits.append(file+'.'+dict_dataset[dataset]+'SVM.'+str(scr+1)+'.'+str(split+1))
                
                #make predictions for each set with the corresponding parameters, only tranied with true labels scr=0
                if scr ==0:
                    clf.fit(X,Y)
                    pred = clf.predict(X_ambi)
                    pred_proba = clf.predict_proba(X_ambi)[:,1]
                    count = [1 for c in range(X_ambi.shape[0])]
                    df_curr_pred = pd.DataFrame()
                    df_curr_pred['drug'] = ambi_list
                    df_curr_pred['count'] = count
                    df_curr_pred['predicted_class'] = pred
                    df_curr_pred['predicted_prob'] = pred_proba
                    #aggregate predictions of same drug
                    df_curr_pred_agg = df_curr_pred.groupby(df_curr_pred['drug']).sum()
                    df_curr_pred_agg.reset_index(inplace=True)
                    
                    df_pred = df_pred_temp.copy()
                    for num,drug in enumerate(df_pred['Compound.Name'].tolist()):
                        count = df_curr_pred_agg[df_curr_pred_agg['drug']==drug]['count'].item()
                        #1st criterion majority of predictions
                        pred_count = df_curr_pred_agg[df_curr_pred_agg['drug']==drug]['predicted_class'].item()
                        #2nd criterion tendency of predicted probabilities
                        proba_count = df_curr_pred_agg[df_curr_pred_agg['drug']==drug]['predicted_prob'].item()
                        if (pred_count/count) > 0.5:
                            df_pred.iloc[num,1] == 'DILI'
                        elif (pred_count/count) < 0.5:
                            df_pred.iloc[num,1] == 'NoDILI'
                        elif (pred_count/count) == 0.5:
                            if (proba_count/count) > 0.5:
                                df_pred.iloc[num,1] == 'DILI'
                            elif (proba_count/count) < 0.5:
                                df_pred.iloc[num,1] == 'NoDILI'
                    
                    predictions.append(df_pred)
                    predictions_ident.append(file+'.'+dict_dataset[dataset]+'RF.'+str(i+1)+'.'+str(j+1))
                    
                    clf2.fit(X,Y)
                    pred = clf2.predict(X_ambi)
                    pred_proba = clf2.predict_proba(X_ambi)[:,1]
                    count = [1 for c in range(X_ambi.shape[0])]
                    df_curr_pred = pd.DataFrame()
                    df_curr_pred['drug'] = ambi_list
                    df_curr_pred['count'] = count
                    df_curr_pred['predicted_class'] = pred
                    df_curr_pred['predicted_prob'] = pred_proba
                    #aggregate predictions of same drug
                    df_curr_pred_agg = df_curr_pred.groupby(df_curr_pred['drug']).sum()
                    df_curr_pred_agg.reset_index(inplace=True)
                    
                    df_pred = df_pred_temp.copy()
                    for num,drug in enumerate(df_pred['Compound.Name'].tolist()):
                        count = df_curr_pred_agg[df_curr_pred_agg['drug']==drug]['count'].item()
                        pred_count = df_curr_pred_agg[df_curr_pred_agg['drug']==drug]['predicted_class'].item()
                        proba_count = df_curr_pred_agg[df_curr_pred_agg['drug']==drug]['predicted_prob'].item()
                        if (pred_count/count) > 0.5:
                            df_pred.iloc[num,1] == 'DILI'
                        elif (pred_count/count) < 0.5:
                            df_pred.iloc[num,1] == 'NoDILI'
                        elif (pred_count/count) == 0.5:
                            if (proba_count/count) > 0.5:
                                df_pred.iloc[num,1] == 'DILI'
                            elif (proba_count/count) < 0.5:
                                df_pred.iloc[num,1] == 'NoDILI'
                    
                    predictions.append(df_pred)
                    predictions_ident.append(file+'.'+dict_dataset[dataset]+'SVM.'+str(i+1)+'.'+str(j+1))
 


        
#export predictions in zip file, code snippet only required to make predictions for FDA test set
"""
zf = zipfile.ZipFile('Predictions_GEX.zip', mode='w')
for i,j in zip(predictions,predictions_ident):
    name = j.replace('.','_')
    i.to_csv(path_or_buf='Predictions_GEX/'+name+'.txt',sep=',',index=False)
    zf.write('Predictions_GEX/'+name+'.txt')
    #os.remove('/ECFP/Predictions_GEX/'+name+'.txt')
zf.close()
"""

#export test performances
df_ts = pd.DataFrame()
df_ts['splits'] = ts_splits
df_ts['accuracy'] = acc_ts
df_ts['balanced_accuracy'] = bal_acc_ts
df_ts['recall'] = rec_ts
df_ts['precision'] = pre_ts
df_ts['AU_Prec_Rec_Curve'] = aupr_ts
df_ts['ROC_AUC'] = rocauc_ts
df_ts['MCC'] = mcc_ts

df_ts.to_csv('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameter/MCNC_ts_scores_GEX_yscr.csv',sep=',',index=False)


#export cv performances
df_cv = pd.DataFrame()
df_cv['splits'] = cv_splits
df_cv['accuracy'] = acc
df_cv['balanced_accuracy'] = bal_acc
df_cv['recall'] = rec
df_cv['precision'] = pre
df_cv['AU_Prec_Rec_Curve'] = aupr
df_cv['ROC_AUC'] = rocauc
df_cv['MCC'] = mcc

df_cv.to_csv('CAMDA-DILI/Machine_Learning/data/Model_Results_Parameter/MCNC_cv_scores_GEX_yscr.csv',sep=',',index=False)

