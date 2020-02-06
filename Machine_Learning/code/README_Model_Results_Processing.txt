(1)
To re-create the python3.7 environment used in this analysis please see:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments "Building identical conda environments"

Specifically run "conda create --name myenv --file CAMDA_Metrics_conda_env_spec.txt".

(2)
Scripts: Parse_Ambiguous_Results_FDA_Validation_Set.ipynb
Purpose: Converts CAMDA confusion matrices to tabulated .csv file with model metrics for models tested on the FDA validation set. Outputs an "ambiguous_results.csv" per model type e.g. inside each of the model folders

--- ambiguous_results.csv file description per column ---

(i) "model"; this is currently encoded as "dataset_" (all - DILIrank + SIDER, MCNC - DILIrank (-vLessConcern), MCLCNC - DILIrank),  + "algorithm." (ML algorithm used) + "1." (redundant) + "train_test_split." (i.e. count over the outer loop of 10 different train-test splits during model generation).

(ii) all other columns are the metrics of each model when tested
(3)
Scripts: Parse_Ambiguous_Results_CV_and_External_Test.ipynb + Parse_Ambiguous_Results_FDA_Validation_Set.ipynb

Purpose: Calculates mean+-std values for each model type for RF and SVM seperately. Such tables were used to produce the rows in Table S1 of the manuscript


