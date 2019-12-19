#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 10:50:09 2019

@author: mw775
"""

from scipy.cluster.hierarchy import linkage
import pandas as pd

from scipy.spatial.distance import squareform

tani = pd.read_csv('CAMDA-DILI/data/processed_data/Standardization/tanimoto_similarities.txt', delimiter=',', header=None)
tani.columns = ['C1', 'C2', 'Tanimoto']
tani = tani.pivot(index='C1', columns='C2', values='Tanimoto')
#leave ambiguous out
tani = tani.iloc[:920,:920]

#convert similarity matrix in distance matrix
dist_matrix = 1 - tani

#convert a matrix in condensed 1d vector
v = squareform(dist_matrix)

#generate the linkage matrix
Z = linkage(v, method='single')

#df_iter : all built clusters
df_iter = pd.DataFrame(Z)
#create dictionary: cluster --> contained mols
dict_dend = {}
for i in range(920):
    dict_dend[i] = [i]

for i,row in df_iter.iterrows():
    if row[0] <= 919:
        dict_dend[i+920] = [row[0]]
        
    if row[0] > 919:
        dict_dend[i+920] = [j for j in dict_dend[row[0]]]
        
              
    if row[1] <= 919:
        dict_dend[i+920].append(row[1])
        
    if row[1] > 919:
        for j in dict_dend[row[1]]:
            dict_dend[i+920].append(j)
            

#which clusters are present at distance cut-off 0.5?

target_set = set(range(920))
df_Z = pd.DataFrame(Z, columns=['cluster1', 'cluster2', 'distance', 'size'])
for i,row in df_Z.iterrows():
    if row['distance']>= 0.5:
        break
    target_set.remove(int(row['cluster1']))
    target_set.remove(int(row['cluster2']))
    target_set.add(int(i+920))
    
#target set: clusters with distance < .5
#dict_dend: compounds in cluster
#dict_cluster: cluster membership at distance .5
dict_cluster = {}
for cluster in target_set:
    for compound in dict_dend[cluster]:
        dict_cluster[int(compound)] = cluster
        
df_compounds = pd.read_csv('CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv', delimiter=',')
clusters = [dict_cluster[i] for i in range(920)]
df_compounds['cluster'] = clusters
df_compounds.to_csv('CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_excl_ambiguous_cluster.csv', sep=',')
