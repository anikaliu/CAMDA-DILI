This file describes how the challenge data are processed to obtain the files for further analysis. 

The script CAMDA-DILI/code/SIDER_get_inactives_liver.py gets compounds from SIDER with no liver related side effects.

  Inputs: 'CAMDA-DILI/processed_data/SIDER/sider_meddra_all_se.csv', 'CAMDA-DILI/processed_data/SIDER/preferred_terms_liver.txt',
  'CAMDA-DILI/processed_data/SIDER/sider_names.csv'

  Output: 'CAMDA-DILI/processed_data/SIDER/sider_inactives_unfiltered.csv'

The script CAMDA-DILI/code/Standardisation_filtering.py combines the challenge dataset with the compounds from SIDER
and standardises them

  Inputs: 'CAMDA-DILI/processed_data/SIDER/sider_inactives_unfiltered.csv', 'CAMDA-DILI/challenge_data/DILIrank_1.csv',      
          '1136compounds.csv'
  Outputs: 'CAMDA-DILI/processed_data/Standardization/standardized_compounds_incl_ambiguous.csv',
           'CAMDA-DILI/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv'
