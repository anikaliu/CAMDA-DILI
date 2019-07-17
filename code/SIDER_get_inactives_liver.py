#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 12:12:57 2019

@author: mw775
"""
#this script generates a dataframe containing all drugs of SIDER with no reported liver side effects

import pandas as pd
import re
import requests
import json
import time
from rdkit import Chem

col = ['stitch_id_flat', 'stitch_id_stereo', 'umls_id_label', 'meddra_conc_type', 'umls_id_medra','side_effect']

#import file containing all drug - side effect relations in SIDER (23/05/19)
df = pd.read_csv('CAMDA-DILI/processed_data/SIDER/sider_meddra_all_se.csv', delimiter='\t', header=None, names=col)


#map STITCH to Pubchem

p = re.compile('CID0*([0-9]+)')



def to_pubchem(stitch_id):
        
    compound_id = p.match(stitch_id).group(1)
      
      
    pubchem_dict['CID' + compound_id] = {'original_id': stitch_id}
    
    




#import file containing liver related preferred terms (pt)

a = []
f = open('CAMDA-DILI/processed_data/SIDER/preferred_terms_liver.txt', 'r')
f1 = f.readlines()
for i in f1:
    if i != '' and i!='\n' and i!='\t\n':
        #print(i)
        a.append(i)
f.close

#list of all luver related pts
dili_pt = [i[1:-1] for i in a]

#get ids with dili se
dili_ids = []
dili_ids = df[df['side_effect'].isin(dili_pt)]['stitch_id_stereo'].unique().tolist()

#dataframe containing only drugs with no side effects contained in dili_pt
no_dili_df = df[~df['stitch_id_stereo'].isin(dili_ids)]

pubchem_dict = dict()
for compound in no_dili_df['stitch_id_stereo'].unique():
    to_pubchem(compound)

#remove CID
pcid = [i[3:] for i in pubchem_dict]


#create df stitch_stereo - pubchem
stitch_stereo = [pubchem_dict[i]['original_id'] for i in pubchem_dict]
df3 = pd.DataFrame()
df3['stitch_stereo'] = stitch_stereo
df3['pubchem'] = pcid


#get Inchi keys

base = 'https://www.ebi.ac.uk/unichem/rest/'

molecule_dict = dict()
errors = []

def find_inchi(compound_id):
    """For a PubChem compound id, retrieve the chembl_id via Unichem
    compound_id -- str, PubChem CID"""
    
    response = requests.get(base + 'structure/{}/22'.format(compound_id)) # source 22 is PubChem
    
    try:
        assert response.status_code == 200
    except AssertionError:
        errors.append(compound_id)
        return
    
    result = json.loads(response.content.decode())
    
    try:
        inchi_key = result[0]['standardinchikey']
        inchi = result[0]['standardinchi']
    except IndexError:
        errors.append(compound_id)
        return
    time.sleep(0.2)
    
    molecule_dict[compound_id] = {'original_id': compound_id, 'inchi_key': inchi_key, 'inchi' : inchi}
    
for i in pcid:
    find_inchi(i)


#dataframe sider_name, pubchem id, inchi key, smiles
df_ni = pd.DataFrame()
df_ni['pubchem_cid'] = [i for i in molecule_dict]
df_ni['inchi_key'] = [molecule_dict[i]['inchi_key'] for i in molecule_dict]
df_ni['smiles'] = [Chem.MolToSmiles(Chem.inchi.MolFromInchi(molecule_dict[i]['inchi'])) for i in molecule_dict]


#get SIDER names
st_st_ni = []
for i in df_ni['pubchem_cid']:
    st_st_ni.append(df3[df3['pubchem'] == i]['stitch_stereo'].iloc[0])
    
st_fl_ni = []
for i in st_st_ni:
    st_fl_ni.append(df[df['stitch_id_stereo'] == i]['stitch_id_flat'].iloc[0])

col2 = ['stitch_id_flat', 'name']
df_names = pd.read_csv('CAMDA-DILI/processed_data/SIDER/sider_names.csv', delimiter='\t', header=None, names=col2)

names = [df_names[df_names['stitch_id_flat'] == i]['name'].iloc[0] for i in st_fl_ni]

#add names to df
df_ni['names'] = names


#export df
df_ni.to_csv(path_or_buf='CAMDA-DILI/processed_data/SIDER/sider_inactives_unfiltered.csv', sep=',')