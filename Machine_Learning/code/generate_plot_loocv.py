#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 10:59:17 2019

@author: moritz
"""

import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.patches as mpatches

#generate file for similarities mcnc, ambiguous
tani_all = pd.read_csv('CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/tanimoto_similarities.txt', delimiter=',', header=None)
tani_all.columns = ['C1', 'C2', 'Tanimoto']
tani_all = tani_all.pivot(index='C1', columns='C2', values='Tanimoto')

tani_mcrows = pd.concat([tani_all.iloc[:174,:174],tani_all.iloc[:174,434:661],tani_all.iloc[:174,920:]],axis=1)
tani_ncrows = pd.concat([tani_all.iloc[434:661,:174],tani_all.iloc[434:661,434:661],tani_all.iloc[434:661,920:]],axis=1)
tani_ambirows = pd.concat([tani_all.iloc[920:,:174],tani_all.iloc[920:,434:661],tani_all.iloc[920:,920:]],axis=1)

tani_mcncambi = pd.concat([tani_mcrows,tani_ncrows,tani_ambirows],axis=0)

#reset indices
d = {}
for i,j in enumerate(tani_mcncambi.columns.values):
    d[j] = i
tani_mcncambi = tani_mcncambi.rename(d,axis=1)
tani_mcncambi = tani_mcncambi.rename(d,axis=0)  

#for each of the ambiguous compute: i) similarity all (0-400) ii) similarity predicted class (0-173/174-400)
sims_all = []
sims_pre_class = []

import numpy as np
#function that gathers similarities for each compound of a group to compounds of another group (inter)
def five_nn_sim_dg(df, idx_compounds_a, idx_compounds_b):
    
    
    #consider only columns of idx_compounds_b for comparison
    df = df[idx_compounds_b]
    
    for i,row in df.iterrows():
        #consider only rows within idx_compounds_a
        if i != idx_compounds_a:
            continue
        l = list(row)
        l.sort(reverse=True)
    return(np.mean(l[0:5]))

labels = []

for i in range(55):
    sims_all.append(five_nn_sim_dg(tani_mcncambi,i+401,list(range(401))))
    
                              
df_ambi_sim = pd.DataFrame()
df_ambi_sim['similarities_all'] = sims_all


df_ambi_dist = pd.DataFrame()
df_ambi_dist['bins'] = ['0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.8-0.9']
df_ambi_dist['counts_sim_all'] = [0 for i in range(10)]




for i,row in df_ambi_sim.iterrows():
    for lim in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
        
        if row['similarities_all'] <= (lim+0.1) and row['similarities_all'] > lim:
            df_ambi_dist.iloc[int(lim*10),1] +=1
            
            
#generate Plot 2a
train_data = pd.read_csv("CAMDA-DILI/Machine_Learning/data/LOOCV_Results/result_loocv_svm_processed.csv")
test_data = df_ambi_dist

train_data['accuracy_sim_all'] = train_data['accuracy_sim_all']*100
train_data['bins'] = ['0.0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1.0']

# Bring some raw data.
frequencies = train_data['accuracy_sim_all'][0:6]
# In my original code I create a series and run on that, 
# so for consistency I create a series from the list.
freq_series = pd.Series.from_array(frequencies)

x_labels = train_data['bins'][0:6]

# Plot the figure.
plt.figure(figsize=(7, 10))
ax = freq_series.plot(kind='bar', color='blue', width=0.7) # ,align='edge'
ax.set_xlabel("Mean Tanimoto similarity to 5 nearest neighbours",labelpad=20)
ax.set_ylabel("Correctly classified (%)")
ax.set_xticklabels(x_labels)

rects = ax.patches
plt.rcParams.update({'font.size': 10})
plt.rc('axes', labelsize=20)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('legend', fontsize=20)
# Make some labels.
labels = ["n = %d" % i for i in train_data['counts_sim_all'][0:6]]

for rect, label in zip(rects, labels):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height + 1, label,
            ha='center', va='bottom', size=18)
    
ax.text(0.01, 0.96, '(a)', transform=plt.gca().transAxes, size = 20)
plt.xlim((-0.5,5.5))
plt.xticks(rotation=45)
plt.ylim((0,110))   


plt.savefig("CAMDA-DILI/Machine_Learning/plots/Training_Data_AD.svg")
       

