(1)

Script: CAMDA_CV_and_ET_Result_Parsing.ipynb
Purpose: To parse ECFP4, MD, PT (+ reduced models) performance metrics and output tables containing both cross-validation (CV)
         and external test set (ET) for RF and SVM seperately. Such tables were used to produce the 'ET' rows in Table S1 of                                                                          
         the manuscript

--- processed_data/Models/ file description per column ---

(i) "model"; this is currently encoded as "val." (redundant) + "dataset_" (all - DILIrank + SIDER, MCNC - DILIrank (-vLessConcern), MCLCNC - DILIrank),  + "algorithm_" (ML algorithm used) + "1_" (redundant) + "train_test_split" (i.e. over the outer loop of 10 different train-test splits)

(ii) all other columns are the metrics of each model when tested on the FDA validation set
