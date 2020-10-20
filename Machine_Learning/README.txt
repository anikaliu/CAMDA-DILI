(1)
To re-create the python3.6 environment used in this analysis please see:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file

Specifically run "conda env create -f camda_env.yml"

(A)
Script: models_ECFP.py
Purpose: do hyperparameter optimisation and prediction for the internal test sets (including y-scrambled models) as well as the FDA test set for ECFP models, create splits used for MD and PT models.

(B)
Script: models_MD.py
Purpose: do hyperparameter optimisation and prediction for the internal test sets (including y-scrambled models) as well as the FDA test set for molecular descriptor models.

(C)
Script: models_PT.py
Purpose: do hyperparameter optimisation and prediction for the internal test sets (including y-scrambled models) as well as the FDA test set for protein target models.

(D)
Script: loocv_best_model.py
Purpose: do the Leave-one-out-cross-validation experiment for the best model (SVM, ECFP, MostConcern + NoConcern)

(E)
Script: generate_plots_loocv.py
Purpose: gnereate plots for Figure 2a.
