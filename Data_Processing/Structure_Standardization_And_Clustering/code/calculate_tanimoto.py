import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


#function calculates tanimoto similarity
def tanimoto(fp1, fp2):

    in_1 = 0
    in_2 = 0
    in_both = 0

    for a,b in zip(fp1,fp2):
        if a == 1 and b == 1:
            in_both+=1
        if a == 1:
            in_1+=1
        if b == 1:
            in_2+=1
    return(in_both/(in_1 + in_2 - in_both))

#open file with smiles
df = pd.read_csv('CAMDA-DILI/data/processed_data/Standardization/standardized_compounds_incl_ambiguous.csv', delimiter=',')

#generate list of fps from smiles
smiles = df['standardized_smiles'].tolist()
mols = [Chem.MolFromSmiles(smile) for smile in smiles]
fps_bit = [AllChem.GetMorganFingerprintAsBitVect(mol,2, nBits=2048) for mol in mols]

thefile = open('CAMDA-DILI/data/processed_data/Standardization/tanimoto_similarities.txt', 'w')
#
count_1 = 1
for fp in fps_bit:
    count_2 = 1
    for fp2 in fps_bit:
        line = str(count_1) + ',' + str(count_2) + ',' + str(tanimoto(fp, fp2)) + '\n'
        thefile.write(line)
        count_2+=1
    count_1+=1
	
thefile.close()
