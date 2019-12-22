ML model performances for CV and ET for each model type (ECFP, MD, and PT) are currently in the directory ./{model_type}" labelled as "cv_scores_FP_yscr.csv" and "ts_scores_FP_yscr.csv"

--- "ts_scores_FP_yscr.csv" file description per column ---

(i) "splits"; this is currently encoded as  "dataset." (all - DILIrank + SIDER, MCNC - DILIrank (-vLessConcern), MCLCNC - DILIrank),  + "algorithm." (ML algorithm used) + "true_or_y_scrambled_labels_used." (ranges from 1-11, i.e. 1 standard model built [1], followed by 10 y-scrambled models [2-11]) + "train_test_split." (i.e. over the outer loop of 10 different train-test splits)

(ii) all other columns are the metrics of each model when tested on the external test set (varied per outer train-test split)

--- "cv_scores_FP_yscr.csv" file description per column ---

(i) "splits"; this is currently encoded as  "dataset." (all - DILIrank + SIDER, MCNC - DILIrank (-vLessConcern), MCLCNC - DILIrank),  + "algorithm." (ML algorithm used) + "true_or_y_scrambled_labels_used." (ranges from 1-11, i.e. 1 standard model built [1], followed by 10 y-scrambled models [2-11]) + "train_test_split." (i.e. over the outer loop of 10 different train-test splits) + "cross-validation split" (i.e. inner 5-fold GroupKFold loop)

(ii) all other columns are the metrics of each model when tested on the external test set (varied per outer train-test split)
