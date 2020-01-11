This files describes the training of the ML models and the prediction of the challenge compounds. Also the 5NN analysis is described

The scripts models_ECFP.py, models_gridsearch_MD.py, models_PT.py
contain the training procedure using  as features ECFP, molecular descriptors and predicted protein targets respectively.
The steps are:
-iteration through datasets (DILIrank+SIDER, DILIrank, MCNC)
-correct labels and 10 times y-scrambled
-separation of the set in 10 folds, 10 times CV (=internal test sets)
-for each of the training folds a grid search CV for hyperparameter optimisation, the algorithms are RF and SVM
-the best set of parameters is used for predicting the ambiguous compounds
-save CV scores, best hyperparameters, internal test scores, predictions of ambiguous compounds

Inputs for all scripts:
-'CAMDA-DILI/data/processed_data/Models/ECFP/standardized_compounds_excl_ambiguous_cluster.csv'
-'CAMDA-DILI/data/processed_data/Models/ECFP/standardized_compounds_incl_ambiguous.csv'
-'CAMDA-DILI/processed_data/Models/myname_predictions_no1_TEMPLATE.txt'

Input only MD:
-'CAMDA-DILI/data/processed_data/Models/MD/mol_descriptors_training.csv'
-'CAMDA-DILI/data/processed_data/Models/MD/mol_descriptors_ambiguous.csv'

Input only PT(note the splits that are imported are created in teh ECFP script, this is necessary since the PT dataset is slightly
smaller than the whole set):
-'CAMDA-DILI/data/processed_data/Models/PT/training_set_predicted_targets.txt'
-'CAMDA-DILI/processed_data/Models/PT/ambiguous_predicted_targets.txt'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/all_ttsplits_tr.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/all_ttsplits_te.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/mclcnc_ttsplits_tr.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/mclcnc_ttsplits_te.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/mcnc_ttsplits_tr.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/mcnc_ttsplits_te.pkl'


Output ECFP:
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/all_ttsplits_tr.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/all_ttsplits_te.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/mclcnc_ttsplits_tr.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/mclcnc_ttsplits_te.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/mcnc_ttsplits_tr.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/splits/mcnc_ttsplits_te.pkl'
-'CAMDA-DILI/data/processed_data/Models/ECFP/Predictions_ECFP.zip'
-'CAMDA-DILI/data/processed_data/Models/ECFP/test_scores_ECFP.csv'
-'CAMDA-DILI/data/processed_data/Models/ECFP/cv_scores_ECFP.csv'
-'CAMDA-DILI/data/processed_data/Models/ECFP/best_parameters_ECFP.pkl'

Output MD:
-'CAMDA-DILI/data/processed_data/Models/MD/Predictions_MD.zip'
-'CAMDA-DILI/data/processed_data/Models/MD/test_scores_MD.csv'
-'CAMDA-DILI/data/processed_data/Models/MD/cv_scores_MD.csv'
-'CAMDA-DILI/data/processed_data/Models/MD/best_parameters_MD.pkl'

Output PT:
-'CAMDA-DILI/data/processed_data/Models/PT/Predictions_PT.zip'
-'CAMDA-DILI/data/processed_data/Models/PT/test_scores_PT.csv'
-'CAMDA-DILI/data/processed_data/Models/PT/cv_scores_PT.csv'
-'CAMDA-DILI/data/processed_data/Models/PT/best_parameters_PT.pkl'
-'CAMDA-DILI/data/processed_data/Models/PT/svm_feature_coefficients.csv'
-'CAMDA-DILI/processed_data/Models/PT/rf_feature_importance.csv'


The script  CAMDA-DILI/code/loocv_best_model.py takes the compounds file 
('CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv') as input and creates the file 
CAMDA-DILI/data/processed_data/LOOCV/result_loocv_svm_processed.csv which contains the results of the LOOCV experiment.

The script  CAMDA-DILI/code/ML_training/generate_plots_loocv.py takes as input the results of the LOOCV experiment
(CAMDA-DILI/data/processed_data/LOOCV/result_loocv_svm_processed.csv) and the Tanimoto similarities
(CAMDA-DILI/data/processed_data/Standardization/tanimoto_similarities.txt) to generate the plots 2a 
(CAMDA-DILI/data/processed_data/Plots/Training_Data_AD.svg) and 2b (CAMDA-DILI/dataprocessed_data/Plots/Training_vs_validation_AD.svg) 
