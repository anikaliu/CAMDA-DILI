(1)
To re-create the python3.6 environment used in this analysis please see:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file

Specifically run "conda env create -f camda_env.yml"

(A)
Script: SIDER_get_inactives_liver.py
Purpose: get inactives compounds (no liver side effects reported) from SIDER database

(B)
Script: Standardisation_filtering.py
Purpose: combine DILIrank and SIDER compounds, standardise chemical structures and filter out compounds (inorganics, polymers, mixtures, out of weight range, duplicates)

(C)
Script: calculate_tanimoto.py
Purpose: calculate Tanimoto similarities between all compounds for clustering

(D)
Script: Hierarchical_Clustering.py
Purpose: do hierarchical clustering based on Tanimoto similarities

(E)
Script: Generation_md.py
Purpose: calculate molecular descriptors for compounds using the Mordred package

(F)
For this step installation of Pidgin Version 3 (https://pidginv3.readthedocs.io/en/latest/) is required. The file running_pidgin.txt describes the used input parameters.
Purpose: predict protein targets for the compounds
