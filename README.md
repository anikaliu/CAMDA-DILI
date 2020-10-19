# CAMDA 2019 Drug Safety Challenge

Title: [Prediction and mechanistic analysis of Drug-Induced Liver Injury (DILI) based on chemical structure](https://www.researchsquare.com/article/rs-16599/v1)

Authors: *Anika Liu, Moritz Walter, Peter Wright, Aleksandra Maria Bartosik, Daniela Dolciami, Abdurrahman Elbasir, Hongbin Yang, Andreas Bender*

Here, we provide supporting code and instructions to reproduce the models and results presented in our manuscript. The methods sub-sections of the manuscript are linked to the corresponding repository subfolders below:

* **Data Preparation**: All files and instructions to filter compounds, compute molecular descriptors and cluster compounds based on Tanimoto similarity of ECFP4 FP can be found in the GitHub repository subfolder 'Data_Processing'
* **Model Generation**: All files required to reproduce the results of this section can be found in the GitHub repository subfolder 'Machine_Learning'. Instructions are provided for [model generation](https://github.com/anikaliu/CAMDA-DILI/blob/master/Machine_Learning/code/README_Model_Generation.txt), and [processing of model results](https://github.com/anikaliu/CAMDA-DILI/blob/master/Machine_Learning/code/README_Model_Results_Processing.txt).
* **Interpretation of Protein Targets**: All files required to reproduce the results of this section can be found in the GitHub repository subfolder 'Machine_Learning' including the corresponding [R Markdown file](https://github.com/anikaliu/CAMDA-DILI/Machine_Learning/code/proteintarget_enrichment.nb.html).
* **Structural Alerts**: All files and instructions to compute and analyse substructures can be found in the GitHub repository subfolder 'Substructural Analysis'
