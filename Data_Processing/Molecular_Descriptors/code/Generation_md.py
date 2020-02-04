#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 10:46:21 2019

@author: mw775
"""
#generate molecular descriptors for training set and ambiguous compounds

import pandas as pd
import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors

df = pd.read_csv('CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv', delimiter=',', index_col=0)


smis = df['standardised_smiles'].tolist()
mols = [Chem.MolFromSmiles(smile) for smile in smis]



calc = Calculator(descriptors, ignore_3D=True)
df_descr = calc.pandas(mols)


df_all = pd.read_csv('CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_incl_ambiguous.csv', delimiter=',')


df_ambis = pd.read_csv('CAMDA-DILI/data/challenge_data/myname_predictions_no1_TEMPLATE.txt'', delimiter=',')
ambis = df_ambis['Compound.Name'].tolist()
ambis_indices = []
for i, row in df_all.iterrows():
    if row['Compound Name'] in ambis and row['vDILIConcern'] == 'Ambiguous DILI-concern':
        ambis_indices.append(i)

df_ambis = df_all.iloc[ambis_indices,:]

smis_a = df_ambis['standardised_smiles'].tolist()
mols_a = [Chem.MolFromSmiles(smile) for smile in smis_a]
df_descr_a = calc.pandas(mols_a)

objs = [df_descr,df_descr_a]
df_d = pd.concat(objs)

#keep only descriptors that could be calculated for all drugs
df_d= df_d.replace([np.inf, -np.inf], np.nan)
df_d = df_d.dropna(axis='columns')
df_d = df_d.select_dtypes(include=['float64', 'int64', 'bool'])
df_d = df_d*1
df_d = df_d.astype('float64')

#export descriptors
df_exp1 = df_d.iloc[:923,:]
df_exp1.to_csv(path_or_buf='CAMDA-DILI/data/processed_data/Molecular_Descriptors/mol_descriptors_training.csv', sep=',')

df_exp2 = df_d.iloc[923:,:]
df_exp2.to_csv(path_or_buf='CAMDA-DILI/data/processed_data/Molecular_Descriptors/mol_descriptors_ambiguous.csv', sep=',')
