(1)
To re-create the python2.7 environment used in this analysis please see:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments "Building identical conda environments"

Specifically run "conda create --name myenv --file CAMDA_SAs_conda_env_spec.txt".

(2)
File order - base file is "compounds_added_inactives.csv":

(A) 
Script: Prepare_DILI_Dataset_for_SARpy_and_MoSS.ipynb
Purpose: Generating DILIrank(-vLessConcern) dataset

(Bi) 
Script: 
Purpose: Generating MoSSpy SAs

(Bii) 
Script: DILI_SUBSTRUCTURES_USE_ALL.ipynb
Purpose: Generating SARpy SAs, then concatenating with MoSS and literature (Liu et al. (2015)) to generate SA metrics e.g. precision and % coverage in DILI compounds

(C)
Script: DrugBank_Approved_Analysis.ipynb
Purpose: Check prescence of all SA in DrugBank v 5.1.4 compounds (standardised in-script)

