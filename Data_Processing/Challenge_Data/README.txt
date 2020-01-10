This file describes how the challenge data are processed to obtain the files for further analysis. Also it describes how molecular descriptors are calculated and protein targets are predicted.

The script CAMDA-DILI/Data_Processing/SIDER_Retrieval/code/SIDER_get_inactives_liver.py gets compounds from SIDER with no liver related side effects.

  Inputs: 'CAMDA-DILI/Data_Processing/SIDER_Retrieval/data/sider_meddra_all_se.csv',
          'CAMDA-DILI/Data_Process/SIDER_Retrieval/data/preferred_terms_liver.txt',
          'CAMDA-DILI/Data_Processing/SIDER_Retrieval/data/sider_names.csv'

  Output: 'CAMDA-DILI/Data_Processing/SIDER_Retrieval/data/sider_inactives_unfiltered.csv'

The script CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/code/Standardisation_filtering.py combines the challenge dataset with the compounds from SIDER and standardises them.

  Inputs: 'CAMDA-DILI/Data_Processing/SIDER_Retrieval/data/sider_inactives_unfiltered.csv',
          'CAMDA-DILI/Data_Processing/challenge_data/DILIrank_1.csv',
          'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/filtered_compounds.csv'
  
  Outputs: 'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_incl_ambiguous.csv',
           'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_excl_ambiguous.csv'
           
The script CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/code/calculate_tanimoto.py calculates Tanimoto similarities between all compounds.
  
  Input: 'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_incl_ambiguous.csv'
  
  Output: 'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/tanimoto_similarities.txt
       
The script CAMDA-DILI/Data_Processing/Molecular_Descriptors/code/Generation_md.py calculates Mordred descriptors for all compounds
  
  Inputs: 'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_excl_ambiguous.csv'
          'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/standardized_compounds_incl_ambiguous.csv'
          'CAMDA-DILI/Data_Processing/challenge_data/myname_predictions_no1_TEMPLATE.txt'
  
  Outputs: 'CAMDA-DILI/Data_Processing/Molecular_Descriptors/data/mol_descriptors_training.csv'
           'CAMDA-DILI/Data_Processing/Molecular_Descriptors/data//mol_descriptors_ambiguous.csv'
           
The script CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/code/Hierarchical_Clustering.py performs hierarchical Clustering for the training set (without ambiguous) based on Tanimoto similarities.
  
  Inputs:'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering/data/tanimoto_similarities.txt
           'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering   /data/standardized_compounds_excl_ambiguous.csv'
          
  Output: 'CAMDA-DILI/Data_Processing/Structure_Standardization_And_Clustering   /data/standardized_compounds_excl_ambiguous_cluster.csv'
  
Target prediction using PIDGIN is done according to the file  CAMDA-DILI/Data_Processing/Target_Prediction/code/running_pidgin.txt 
  
  Inputs: 'CAMDA-DILI/Data_Processing/Target_Prediction/data/training_set_pidgin.smi'
          'CAMDA-DILI/Data_Processing/Target_Prediction/data/ambiguous_pidgin.smi'
        
  Outputs: 'CAMDA-DILI/Data_Processing/Target_Prediction/data/training_set_predicted_targets.txt'
           'CAMDA-DILI/Data_Processing/Target_Prediction/data/ambiguous_predicted_targets.txt '
