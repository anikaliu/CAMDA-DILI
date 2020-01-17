#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:56:42 2020

@author: moritz
"""
#check which compounds are removed by the standardization
from __future__  import print_function, division, absolute_import
from rdkit.Chem.Descriptors import ExactMolWt

import warnings


import pandas as pd
 
from rdkit import Chem


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


input_dili = pd.read_csv('CAMDA-DILI/Data_Processing/challenge_data/DILIrank_1.csv')

input_sider = pd.read_csv('CAMDA-DILI/Data_Processing/SIDER_Retrieval/data/sider_inactives_unfiltered.csv')


#which PubChem_CIDs from the inputs are not in output
cid_input_dili = set(input_dili['PubChem_CID'].tolist())
cid_input_sider = set(input_sider['pubchem_cid'].tolist())

#combine inputs

#remove entries from sider which are already in dili
olap_dili_sider = cid_input_dili.intersection(cid_input_sider)

#get indices of overlap
indices_to_drop = [] 
for i,row in input_sider.iterrows():
    if row['pubchem_cid'] in olap_dili_sider:
        indices_to_drop.append(i)
        
input_sider.drop(labels=indices_to_drop,inplace=True)

input_sider['vDILIConcern'] = ['sider_inactive' for i in range(448)]

#match columns of sider
input_sider.rename(mapper={'pubchem_cid':'PubChem_CID','names':'Compound Name','smiles':'SMILES'}, axis=1, inplace=True)

joined_inputs = pd.concat(objs=[input_dili,input_sider],join='inner')
joined_inputs.reset_index(drop=True,inplace=True)



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
        


#remove compounds that can't form a rdkit mol
indices_to_drop = []
smiles1 = joined_inputs['SMILES'].tolist()
for i,smi in enumerate(smiles1):
    m = Chem.MolFromSmiles(smi)
    if m == None:
        indices_to_drop.append(i)
        
joined_inputs.drop(indices_to_drop,inplace=True)

smiles1 = joined_inputs['SMILES'].tolist()
can_smis1 = [Chem.MolToSmiles(Chem.MolFromSmiles(smi),True) for smi in smiles1]
joined_inputs['can_smis'] = can_smis1 

#drop mols with MW>1500


mol_wts = []
for smi in can_smis1:
    m = Chem.MolFromSmiles(smi)
    mol_wts.append(ExactMolWt(m))

joined_inputs['mol_weight'] = mol_wts
inputs = joined_inputs[joined_inputs['mol_weight']<1500]
inputs.reset_index(drop=True,inplace=True)

input_smis = inputs['can_smis'].tolist()
input_mols = [Chem.MolFromSmiles(smi) for smi in input_smis]
output_smis = []
for m in input_mols:
    output_smis.append(standard(m))
    

#the following compounds are filtered out (indices of input df)

#mixtures: 5,90,186,461,512,549,671,814,815,966,999,1011,1154,1165,1169,1172,1174,1191,1197,1215,1218,1234,1246,1254,1257,1262,1265,1273,1281,1286,1287,1302,1321,1328,1351,1360,1363

#outside of weight range (>1000 Da): 114,115,152,227,342,418,448,465,468,497,562,606,623,673,684,710,714,718,720,799,820,849,854,857,862,927,1016,1065,1114,1139,1140,1143,1204,1214,1247,1248,1252,1255,1256,1305,1311,1331,1353,1359,1370,

#removed by Flatkinson Standardiser: 106,156,483,532,730,933,952,954,957,959,961,962,1242,944

#inorganic/organometallic: 161,192,462,496,537,570,614,626,672,674,677,792,846,848,947,949,950,953,965,967,968,971,972,973,974,975,1075,1112,1142,1152,1155,1164,1166,1168,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1193,1194,1195,1198,1217,1241,1258,1259,1264,1272,1278,1284,1306,1315,1316,1317,1322,1323,1342,1345,1347,1352,1354,1362

#duplicates or different enantiomers: 308,500,501,503,574,879,945,946,951,955,956,970,978,984,994,1000,1003,1018,1020,1024,1035,1048,1049,1052,1055,1059,1060,1061,1066,1076,1083,1085,1088,1091,1093,1097,1098,1111,1124,1126,1127,1128,1130,1131,1137,1151,1158,1163,1192,1199,1201,1206,1208,1222,1226,1276,1277,1290,1294,1301,1313

#polymers: 332,644,645,1104,1283

indices_to_remove = [5,90,186,461,512,549,671,814,815,966,999,1011,1154,1165,1169,1172,1174,1191,1197,1215,1218,1234,1246,1254,1257,1262,1265,1273,1281,1286,1287,1302,1321,1328,1351,1360,1363,114,115,152,227,342,418,448,465,468,497,562,606,623,673,684,710,714,718,720,799,820,849,854,857,862,927,1016,1065,1114,1139,1140,1143,1204,1214,1247,1248,1252,1255,1256,1305,1311,1331,1353,1359,1370,944,106,156,483,532,730,933,954,957,959,961,962,1242,161,192,462,496,537,570,614,626,672,674,677,792,846,848,947,949,950,953,965,967,968,971,972,973,974,975,1075,1112,1142,1152,1155,1164,1166,1168,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1193,1194,1195,1198,1217,1241,1258,1259,1264,1272,1278,1284,1306,1315,1316,1317,1322,1323,1342,1345,1347,1352,1354,1362,308,500,501,503,574,879,945,946,951,955,956,970,978,984,994,1000,1003,1018,1020,1024,1035,1048,1049,1052,1055,1059,1060,1061,1066,1076,1083,1085,1088,1091,1093,1097,1098,1111,1124,1126,1127,1128,1130,1131,1137,1151,1158,1163,1192,1199,1201,1206,1208,1222,1226,1276,1277,1290,1294,1301,1313,332,644,645,1104,1283,952]

#delete compounds from standardised output
for i in sorted(indices_to_remove, reverse=True):
    del output_smis[i]
    
#delete compounds from input_df
inputs.drop(labels=indices_to_remove,inplace=True)


#keep first fragment
for en,i in enumerate(output_smis):
    if type(i)==list:
        output_smis[en] = i[0]
        
inputs['standardised_smiles'] = output_smis
inputs.drop(labels=['can_smis','mol_weight'],axis=1,inplace=True)
inputs.drop(labels=['SMILES'],axis=1,inplace=True)
cols = ['PubChem_CID', 'Compound Name', 'standardised_smiles', 'vDILIConcern']
inputs = inputs[cols]

#export file
inputs.to_csv('CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_incl_ambiguous.csv', sep=',',index=False)

#file without ambiguous
df_nonambi =inputs[inputs['vDILIConcern'] != 'Ambiguous DILI-concern']
df_nonambi.to_csv(path_or_buf='CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_excl_ambiguous.csv', sep=',',index=False)
