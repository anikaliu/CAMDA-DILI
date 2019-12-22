(1)
To re-create the python3.7 environment used in this analysis please see:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments "Building identical conda environments"

Specifically run "conda create --name myenv --file CAMDA_Metrics_conda_env_spec.txt".

(2)
Script: Result_Parsing.ipynb
Purpose: Converts confusion matrices to tabulated .csv file with model metrics (all models within nested CV loop). Outputs an "ambiguous_results.csv" per model type e.g. inside each of the model folders: ECFP4, ECFP4_red etc.

--- ambiguous_results.csv file description per column ---

(i) "model"; this is currently encoded as "val." (redundant) + "dataset_" (all - DILIrank + SIDER, MCNC - DILIrank (-vLessConcern), MCLCNC - DILIrank),  + "algorithm_" (ML algorithm used) + "true_or_y_scrambled_labels_used." (ranges from 1-11, i.e. 1 standard model built [1], followed by 10 y-scrambled models [2-11]) + "train_test_split." (i.e. over the outer loop of 10 different train-test splits)

(ii) all other columns are the metrics of each model when tested on the FDA validation set

(3)
Script: Parse_Ambiguous_Results.ipynb
Purpose: Calculates mean+-std values for each model type for the FDA validation set for RF and SVM seperately. Such tables were     used to produce the 'FDA' rows in Table S1 of the manuscript


