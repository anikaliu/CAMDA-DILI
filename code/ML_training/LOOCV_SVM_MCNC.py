# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from sklearn.svm import SVC

#import data table
df = pd.read_csv('CAMDA-DILI/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv', delimiter=',')

#generate FPs
smis = df['standardized_smiles'].tolist()
mols = [Chem.MolFromSmiles(smile) for smile in smis]
fps_bit = [Chem.GetMorganFingerprintAsBitVect(mol,2, nBits=2048) for mol in mols]

#store FP bits as input in X
X = np.empty([920,2048])
for i,j in enumerate(fps_bit):
    for k,l in enumerate(list(j)):
        X[i,k] = l
X = np.concatenate((X[:174,:],X[434:661,:]))
       
# labels (0,1) in Y
Y = np.empty([920,1])
for i,j in enumerate(df['vDILIConcern']):
    if j == 'vMost-DILI-Concern' or j == 'vLess-DILI-Concern':
        Y[i,0] = 1
        
    if j == 'vNo-DILI-Concern' or j == 'sider_inactive':
        Y[i,0] = 0
Y = np.squeeze(Y)
Y= np.concatenate((Y[:174],Y[434:661]))

#import similarities
#tani = pd.read_csv('CAMDA-DILI/processed_data/Standardization/tanimoto_similarities.txt', delimiter=',', header=None)
tani = pd.read_csv('tanimoto_similarities.txt', delimiter=',', header=None)
tani.columns = ['C1', 'C2', 'Tanimoto']
tani = tani.pivot(index='C1', columns='C2', values='Tanimoto')

tani_mcnc = pd.concat([tani.iloc[:174,:174],tani.iloc[:174,434:661]],axis=1)
tani_mcnc = pd.concat([tani_mcnc,pd.concat([tani.iloc[434:661,:174],tani.iloc[434:661,434:661]],axis=1)],axis=0)

d = {}
for i,j in enumerate(tani_mcnc.columns.values):
    d[j] = i

tani_mcnc = tani_mcnc.rename(d,axis=1)
tani_mcnc = tani_mcnc.rename(d,axis=0)    

#function that gathers similarities for each compound of a group to compounds of another group (inter)

def five_nn_sim_dg(df, idx_compounds_a, idx_compounds_b):
    sims = []
    
    #consider only columns of idx_compounds_b for comparison
    df = df[idx_compounds_b]
    
    for i,row in df.iterrows():
        #consider only rows within idx_compounds_a
        if i != idx_compounds_a:
            continue
        l = list(row)
        l.sort(reverse=True)
        sims.append(np.mean(l[0:5]))
    return(sims)

#function to check how many of the 5nn in positives
def five_nn_in_pos(df,idx,idx_to_compare,idx_pos):
    l1 = df[idx].tolist()
    d = {}
    #dictionary(index in table : similarity)
    for i,j in enumerate(l1):
        if i == idx:
            continue
        d[i] = j
    sims_tuples = sorted(d.items(),reverse=True, key=lambda x: x[1])
    indices = [elem[0] for elem in sims_tuples]
    counter = 0
    for i in range(5):
        if indices[i] in list(range(174)):
            counter+=1
    return(counter)
    
    
sim_all = []
sim_1 = []
sim_0 = []
nn_1 = []
nn_0 = []

true_label = []
pr_label = []
correct = []
probability = []

clf = SVC(kernel='linear', C=0.1,class_weight='balanced', random_state=22,probability=True)


for i in range(X.shape[0]):
    sim_all.append(five_nn_sim_dg(tani_mcnc,i,set(range(401))-set([i]))[0])
    sim_1.append(five_nn_sim_dg(tani_mcnc,i,set(range(174))-set([i]))[0])
    sim_0.append(five_nn_sim_dg(tani_mcnc,i,set(range(174,401))-set([i]))[0])
    nn_1.append(five_nn_in_pos(tani_mcnc,i,set(range(401))-set([i]),range(174)))
    nn_0.append(5-nn_1[i])
    
    Xi = X[list(set(range(401))-set([i])),:]
    Yi = Y[list(set(range(401))-set([i]))]
    clf.fit(Xi,Yi)
    
    true_label.append(Y[i])
    pr_label.append(clf.predict(X[i,:].reshape(1,-1))[0])
    correct.append(((true_label[i]==pr_label[i])*1))
    probability.append(clf.predict_proba(X[i,:].reshape(1,-1))[0])
    
    print(i,'/401')
    
df_results = pd.DataFrame()
df_results['similarities_all'] = sim_all
df_results['similarities_DILI'] = sim_1
df_results['similarities_noDILI'] = sim_0
df_results['nn_DILI'] = nn_1
df_results['nn_noDILI'] = nn_0
df_results['true_label'] = true_label
df_results['predicted_label'] = pr_label
df_results['correct'] = correct
df_results['probability'] = probability



df = df_results.copy()

#sim true and false
df_result = pd.DataFrame()
df_result['bins'] = ['0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.8-0.9']
df_result['counts_sim_all'] = [0 for i in range(10)]
df_result['correct_sim_all'] = [0 for i in range(10)]
df_result['accuracy_sim_all'] = [None for i in range(10)]
df_result['counts_sim_true_class'] = [0 for i in range(10)]
df_result['correct_sim_true_class'] = [0 for i in range(10)]
df_result['accuracy_sim_true_class'] = [None for i in range(10)]
df_result['counts_sim_false_class'] = [0 for i in range(10)]
df_result['correct_sim_false_class'] = [0 for i in range(10)]
df_result['accuracy_sim_false_class'] = [None for i in range(10)]


for i,row in df.iterrows():
    for lim in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
        #get bins for whole set
        if row['similarities_all'] <= (lim+0.1) and row['similarities_all'] > lim:
            df_result.iloc[int(lim*10),1] +=1
            if row['correct'] == 1:
                df_result.iloc[int(lim*10),2] +=1
        
        #if true DILI what is similarity to DILI?
        if row['true_label'] == 1:
            if row['similarities_DILI'] <= (lim+0.1) and row['similarities_DILI'] > lim:
                df_result.iloc[int(lim*10),4] +=1
                if row['correct'] == 1:
                    df_result.iloc[int(lim*10),5] +=1
                    
            #...and to NoDILI?
            if row['similarities_noDILI'] <= (lim+0.1) and row['similarities_noDILI'] > lim:
                df_result.iloc[int(lim*10),7] +=1
                if row['correct'] == 1:
                    df_result.iloc[int(lim*10),8] +=1
        
        #if predicted noDILI, what is similarity to NoDILI
        if row['true_label'] == 0:
            if row['similarities_noDILI'] <= (lim+0.1) and row['similarities_noDILI'] > lim:
                df_result.iloc[int(lim*10),4] +=1
                if row['correct'] == 1:
                    df_result.iloc[int(lim*10),5] +=1
            
            #...and to DILI?
            if row['similarities_noDILI'] <= (lim+0.1) and row['similarities_noDILI'] > lim:
                df_result.iloc[int(lim*10),7] +=1
                if row['correct'] == 1:
                    df_result.iloc[int(lim*10),8] +=1        
                    
for i,row in df_result.iterrows():
    if row['counts_sim_all'] > 0:
        df_result.iloc[i,3] = (row['correct_sim_all']/row['counts_sim_all'])
    if row['counts_sim_true_class'] > 0:
        df_result.iloc[i,6] = (row['correct_sim_true_class']/row['counts_sim_true_class'])
    if row['counts_sim_false_class'] > 0:
        df_result.iloc[i,9] = (row['correct_sim_false_class']/row['counts_sim_false_class'])
        
#get for each molecule (0-400) predicted class, correct, sim all, sim true, sim false, sim true - sim false
df_withdif = pd.DataFrame()
df_withdif['predicted_class'] = df['predicted_label']
df_withdif['true_class'] = df['true_label']
df_withdif['correct'] = df['correct']
df_withdif['similarity_all'] = df['similarities_all']
sim_true = []
sim_false = []
for i,row in df.iterrows():
    if row['true_label'] == 1:
        sim_true.append(row['similarities_DILI'])
        sim_false.append(row['similarities_noDILI'])
        
    if row['true_label'] == 0:
        sim_true.append(row['similarities_noDILI'])
        sim_false.append(row['similarities_DILI'])
    
df_withdif['similarity_true'] = sim_true
df_withdif['similarity_false'] = sim_false
df_withdif['true - false'] = [(true - false) for true,false in zip(sim_true,sim_false)]
    
df_result.to_csv('CAMDA-DILI/processed_data/LOOCV/result_loocv_svm_processed.csv', sep=',', index=False)


