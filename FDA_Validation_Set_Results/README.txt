(1)
To re-create the python3.7 environment used in this analysis please see:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments "Building identical conda environments"

Specifically run "conda create --name myenv --file CAMDA_SAs_conda_env_spec.txt".

(2)
Script: Result_Parsing.ipynb
Purpose: Converts confusion matrices to tabulated .csv file with model metrics (all models within nested CV loop). Outputs an "ambiguous_results.csv" per model type e.g. ECFP4, ECFP4_red

(3)
Script: Parse_Ambiguous_Results.ipynb
Purpose: Calculates mean+-std values for each model type for the FDA validation set for RF and SVM seperately. Such tables were     used to produce the 'FDA' rows in Table S1 of the manuscript


