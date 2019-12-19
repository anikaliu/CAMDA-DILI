This file describes how the LOOCV experiment is performed.

The script  CAMDA-DILI/code/LOOCV_SVM_MCNC.py takes the compounds file 
('CAMDA-DILI/processed_data/Standardization/standardized_compounds_excl_ambiguous.csv') as input and creates the file 
CAMDA-DILI/processed_data/LOOCV/result_loocv_svm_processed.csv which contains the results of the LOOCV experiment.

The script  CAMDA-DILI/code/generate_plots_loocv.py takkes as input the results of the LOOCV experiment
(CAMDA-DILI/processed_data/LOOCV/result_loocv_svm_processed.csv) and the Tanimoto similarities
(CAMDA-DILI/processed_data/Standardization/tanimoto_similarities.txt) to generate the plots 2a 
(CAMDA-DILI/processed_data/LOOCV/Training_Data_AD.svg) and 2b (CAMDA-DILI/processed_data/LOOCV/Training_vs_validation_AD.svg)  
