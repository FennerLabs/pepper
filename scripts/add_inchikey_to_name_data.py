import pandas as pd
from pepper_lab.pepper import Pepper
from pepper_lab.datastructurewwtp import DataStructureWWTP
pep = Pepper()
data = DataStructureWWTP(pep)
data.name_data = pd.read_csv('data/wwtp-data/name_data_example_WWTP_combined_data.tsv', sep='\t', index_col=0)
data.name_to_smiles(verbose=True)
data.name_data.to_csv('output/name_data_WWTP_combined_data.tsv', sep='\t', index=False)
