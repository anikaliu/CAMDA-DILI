#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 10:39:28 2019

@author: mw775
"""

from __future__  import print_function, division, absolute_import


import warnings


import pandas as pd
 
from rdkit import Chem
from rdkit.Chem import Descriptors

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = (450, 200)


from standardiser import break_bonds
from standardiser import neutralise
from standardiser import rules
from standardiser import unsalt

break_bonds.logger.setLevel('DEBUG')
neutralise.logger.setLevel('DEBUG')
rules.logger.setLevel('DEBUG')
unsalt.logger.setLevel('DEBUG')

def standard(mol):
    mol = break_bonds.run(mol)
    mol = neutralise.run(mol)
    mol = rules.run(mol)
    mol = neutralise.run(mol)
    #check for fragments
    non_salt_frags = []
    for frag in Chem.GetMolFrags(mol, asMols=True):
        if unsalt.is_nonorganic(frag):
            continue
        if unsalt.is_salt(frag):
            continue
        non_salt_frags.append(Chem.MolToSmiles(frag))
    if len(non_salt_frags) == 1:
        return(non_salt_frags[0])
    else:
        return(non_salt_frags)
        
        
#import DILI rank data
df = pd.read_csv('CAMDA-DILI/Data_Processing/challenge_data/DILIrank_1.csv', delimiter=',')
df = df.drop(df[df.SMILES == '.'].index)
df = df.drop(df[df.SMILES == '0'].index)
smiles1 = df['SMILES'].tolist()
can_smis1 = [Chem.MolToSmiles(Chem.MolFromSmiles(smi),True) for smi in smiles1]
df['can_smis'] = can_smis1        
        
        
input_smis = df['can_smis'].tolist()
input_mols = [Chem.MolFromSmiles(smi) for smi in input_smis]
output_smis = []
for m in input_mols:
    output_smis.append(standard(m))        
        
        


#remove mixtures, inorganic compounds, polymers, salt fragments, duplicates, >1000 Da        
#keep 1st component of the following
keep = [90, 494, 605, 639, 712, 715, 780]
for i,j in enumerate(output_smis):
    if i in keep:
        output_smis[i] = j[0]
df['standardized_smiles'] = output_smis
        






#standardization and filtering for DILI rank complete!

#sider inactives
df2 = pd.read_csv('CAMDA-DILI/Data_Processing/SIDER_Retrieval/data/sider_inactives_unfiltered.csv', delimiter=',')
df2 = df2.drop('Unnamed: 0', axis=1)

df2['MW'] = [Descriptors.ExactMolWt(Chem.MolFromSmiles(smi)) for smi in df2['smiles'].tolist()]

#keep compounds < 1000Da in df_l
df_l = df2.copy()
df_l = df_l[df_l['MW'] < 1000]

input_mols = [Chem.MolFromSmiles(smi) for smi in df_l['smiles'].tolist()]
output_smis = []
for m in input_mols:
    output_smis.append(standard(m))
    
    
#remove bromide, chloride, iodide
for i,j in enumerate(output_smis):
    if type(j) == list and len(j)>1 and j[1]=='[Br-]' :
        output_smis[i] = j[0]
    if type(j) == list and len(j)>1 and j[0]=='[Br-]' :
        output_smis[i] = j[1]
    if type(j) == list and len(j)>1 and j[1]=='[Cl-]' :
        output_smis[i] = j[0]
    if type(j) == list and len(j)>1 and j[0]=='[Cl-]' :
        output_smis[i] = j[1]
    if type(j) == list and len(j)>1 and j[1]=='[I-]' :
        output_smis[i] = j[0]
    if type(j) == list and len(j)>1 and j[0]=='[I-]' :
        output_smis[i] = j[1]
        
#decide which fragment to keep
#remove empty lists
fragments_index = []
for i,j in enumerate(output_smis):
    if type(j) == list and len(j) >1:
        fragments_index.append(i)
        
#125 amphe identic
output_smis[125] = output_smis[125][0]
#322 tranylcypromine racemate, make flat
output_smis[322] = 'NC1CC1c1ccccc1'
#352 identical fragments, keep first
output_smis[352] = output_smis[352][0]
#356 dimenhydrinate, mixture, drop
#397, triphasil, mixture, drop
#400, prussian complex, drop
#429, complex drop
#456, identical fragments, keep first
output_smis[456] = output_smis[456][0]
#465 polymer, drop
#482 mixture, drop
#485 mixture, drop
#497 mixture, drop
#508 mixture, drop
#513 mixture, drop
#516 mixture, drop
#532 mixture, drop
#533 mixture, drop
#551 mixture, drop
#552 drop
#553 mixture, drop
#559 mixture, drop
#581 mixture, drop
#583 complex, drop
#588 complex, drop
#590 pertechnetate, drop
#591 mixture, drop

df_l['standardized_smiles'] = output_smis



df_target = pd.read_csv('CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/Standardization/filtered_compounds.csv', delimiter=',', index_col=0)
target_pc = set(df_target['PubChem_CID'].tolist())
missing_pc = target_pc

for i in df['PubChem_CID']:
    try:
        if int(i) not in target_pc:
            df = df[df['PubChem_CID']!=i]
        if int(i) in target_pc:
            missing_pc = missing_pc - set([int(i)])
    except:
        df = df[df['PubChem_CID']!=i]
        
        
df = df[['PubChem_CID','Compound Name','standardized_smiles','vDILIConcern']]


df_l = df_l[['pubchem_cid','names','standardized_smiles']]
df_l.rename({'pubchem_cid' : 'PubChem_CID', 'names' : 'Compound Name','standardized_smiles' : 'standardized_smiles'},axis=1,inplace=True)
df_l['vDILIConcern'] = 'sider_inactive'

for i in df_l['PubChem_CID']:
    if int(i) not in missing_pc:
        df_l = df_l[df_l['PubChem_CID']!=i]
    if int(i) in missing_pc:
        missing_pc = missing_pc - set([int(i)])

df_final = pd.concat([df,df_l])

df_final.to_csv('CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_incl_ambiguous.csv', sep=',',index=False)


#file without ambiguous
df_nonambi = df_final[df_final['vDILIConcern'] != 'Ambiguous DILI-concern']
df_nonambi.to_csv(path_or_buf='CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_excl_ambiguous.csv', sep=',',index=False)   
