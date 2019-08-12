#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 11:54:27 2019

@author: mw775
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from sklearn.ensemble import RandomForestClassifier


#import data table
df = pd.read_csv('Files/no_ambi_compounds.csv', delimiter=',')
df = df.drop('Unnamed: 0', axis=1)

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
tani = pd.read_csv('Files/Tanimoto/190606results_tanimoto.txt', delimiter=',', header=None)
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

forest = RandomForestClassifier(n_estimators = 750,

                                random_state = 2, max_depth=None,min_samples_leaf=1,

                                class_weight = 'balanced_subsample', bootstrap=True)


for i in range(X.shape[0]):
    sim_all.append(five_nn_sim_dg(tani_mcnc,i,set(range(401))-set([i])))
    sim_1.append(five_nn_sim_dg(tani_mcnc,i,set(range(174))-set([i])))
    sim_0.append(five_nn_sim_dg(tani_mcnc,i,set(range(174,401))-set([i])))
    nn_1.append(five_nn_in_pos(tani_mcnc,i,set(range(401))-set([i]),range(174)))
    nn_0.append(5-nn_1[i])
    
    Xi = X[list(set(range(401))-set([i])),:]
    Yi = Y[list(set(range(401))-set([i]))]
    forest.fit(Xi,Yi)
    
    true_label.append(Y[i])
    pr_label.append(forest.predict(X[i,:].reshape(1,-1)))
    correct.append((true_label[i]==pr_label[i])*1)
    probability.append(forest.predict_proba(X[i,:].reshape(1,-1)))
    
    print(i,'/401', correct[i])
    
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

df_results.to_csv(path_or_buf='Files/result_loocv_mcnc.txt',sep=',')