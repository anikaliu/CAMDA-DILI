(1)
To re-create the python2.7 environment used in this analysis please see:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments "Building identical conda environments"

Specifically run "conda create --name myenv --file CAMDA_SAs_conda_env_spec.txt".

Additionally to replicate the SARpy implementation the following files need to be downloaded:
SARpy.py, SARpy.pyc, SARpytools.py, SARpytools.pyc

(2)
File order - base file is "standardized_compounds_incl_ambiguous.csv" (in CAMDA-DILI/processed_data/Standardization). Additionally the file "DrugBank_5_1_4_approved.csv" is required for step (D):

(A) 
Script: Prepare_DILI_Dataset_for_SARpy_and_MoSS.ipynb
Purpose: (i) Generating DILIrank(-vLessConcern) dataset for SARpy (Most_No_DILIConcern_Dataset_for_SAs.csv is outputted) (ii) Further procecessing to generate .csv file for  MoSS KNIME workflow (moss_input.csv is outputted)

(B)
Script: KNIME workflow in 'KNIME_project MoSS substructures_final.knwf'. Input 'moss_input.csv' via copying data manually into table creator node
Purpose: Generating MoSS SAs (moss_results.csv is outputted)
Note: In "CSV Writer" node; please configure the output location as "GitHub/CAMDA-DILI/Substructural_Analysis/Data/moss_output.csv" directory

(C) 
Script: Generate_SARpy_substructures_and_evaluate_alongside_MoSS_and_Literature
Purpose: Generating SARpy SAs, then concatenating with MoSS and literature (Liu et al. (2015)) to generate SA metrics e.g. precision and % coverage in DILI compounds (Structural_alerts.csv is outputted)

(D)
Script: DrugBank_Approved_Analysis_of_SAs
Purpose: Check prescence of all SA in DrugBank v 5.1.4 compounds (standardised in-script) and concatenate with other metrics  (used for Tables 2, and S4) (Structural_alerts_with_DrugBank.csv is outputted)

(E)
Script: Combined_Plot.R
Purpose: Generate Figure 5
Note: Open in R studio and set working directory to source file location

