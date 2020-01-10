This file describes how the challenge data are processed to obtain the files for further analysis. Also it describes how molecular descriptors are calculated and protein targets are predicted.

The script CAMDA-DILI/code/SIDER_get_inactives_liver.py gets compounds from SIDER with no liver related side effects.

  Inputs: 'CAMDA-DILI/data/processed_data/SIDER/sider_meddra_all_se.csv', 'CAMDA-DILI/data/processed_data/SIDER
  /preferred_terms_liver.txt',
  'CAMDA-DILI/data/processed_data/SIDER/sider_names.csv'

  Output: 'CAMDA-DILI/data/processed_data/SIDER/sider_inactives_unfiltered.csv'

The script CAMDA-DILI/code/Standardisation_filtering.py combines the challenge dataset with the compounds from SIDER
and standardises them

  Inputs: 'CAMDA-DILI/data/processed_data/SIDER/sider_inactives_unfiltered.csv', 'CAMDA-DILI/data/challenge_data
  /DILIrank_1.csv', 'filtered_compounds.csv'
  Outputs: 'CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_incl_ambiguous.csv',
           'CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv'
           
The script calculate_tanimoto.py calculates Tanimoto similarities between all compounds.
  Input: 'CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_incl_ambiguous.csv'
  Output: 'CAMDA-DILI/data/processed_data/Standardization/tanimoto_similarities.txt
       
The script Generation_md.py calculates Mordred descriptors for all compounds
  Inputs: 'CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv'
          'CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_incl_ambiguous.csv'
          'CAMDA-DILI/data/challenge_data/myname_predictions_no1_TEMPLATE.txt'
  
  Outputs: 'CAMDA-DILI/data/processed_data/Molecular_Descriptors/mol_descriptors_training.csv'
           'CAMDA-DILI/data/processed_data/Molecular_Descriptors/mol_descriptors_ambiguous.csv'
           
The script Hierarchical_Clustering.py performs hierarchical Clustering for the training set (without ambiguous) based on Tanimoto similarities.
  Inputs: 'CAMDA-DILI/data/processed_data/Standardization/tanimoto_similarities.txt'
          'CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv'
          
  Output: 'CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_excl_ambiguous_cluster.csv'
  
Target prediction using PIDGIN is done according to the file  CAMDA-DILI/data/processed_data/Target Prediction/running_pidgin.txt 
  
  Inputs: 'CAMDA-DILI/data/processed_data/Target Prediction/training_set_pidgin.smi'
          'CAMDA-DILI/data/processed_data/Target Prediction/ambiguous_pidgin.smi'
        
  Outputs: CAMDA-DILI/data/processed_data/Target Prediction/training_set_predicted_targets.txt
           CAMDA-DILI/data/processed_data/Target Prediction/ambiguous_predicted_targets.txt 
