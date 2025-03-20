# from subprocess import Popen, PIPE
import re
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen
import pubchempy as pcp
import seaborn as sns
import matplotlib.pyplot as plt



class Util:
    @staticmethod
    def name_to_smiles(dataframe):   # todo: can it be non_static, can it be used for any 'data stage'?
        """return a SMILES column in dataframe.
            The names of compounds for which a PubChem entry
            was not found are displayed too
        """
        if 'Name' in dataframe.columns:
            name_list = dataframe.Name
        else:
            name_list = None
            print(dataframe.columns)
            print('Name not in columns ')

        smiles_list = []
        print("names not found in PubChem:")
        for name in name_list:
            smiles = pcp.get_properties(['CanonicalSMILES'], name, 'name', as_dataframe=False)
            if not smiles:
                smiles = np.nan
                print(name)
            else:
                smiles = smiles[0].get('CanonicalSMILES')
            smiles_list.append(smiles)

        # Return a copy of data frame with the new column list
        new_df = dataframe.copy()
        new_df.SMILES = smiles_list

        return new_df

    @staticmethod
    def canonicalize_smiles(smiles, uncharge=False):  # using rdkit
        mol = Chem.MolFromSmiles(smiles)  # creates mol object from SMILES
        if uncharge:
            mol = Chem.MolStandardize.Uncharger().uncharge(mol)  # protonates or deprotonates the mol object
        new_smiles = Chem.rdmolfiles.MolToSmiles(mol)  # converts mol object to canonical SMILES
        can_smiles = Chem.CanonSmiles(new_smiles)
        return can_smiles
    @staticmethod
    def get_clogp(smiles):
        # smiles = str(smiles)
        mol = Chem.MolFromSmiles(smiles)
        try:
            return Crippen.MolLogP(mol)
        except:
            return np.nan
    def remove_stereochemistry(smiles):
        # Convert to mol
        mol = Chem.MolFromSmiles(smiles)
        # Remove stereochemistry
        Chem.rdmolops.RemoveStereochemistry(mol)
        # Convert to back to SMILES
        new_smiles = Chem.MolToSmiles(mol)
        return new_smiles

    @staticmethod
    def smiles_to_inchikey(smiles, uncharge=False):  # using rdkit
        mol = Chem.MolFromSmiles(smiles)  # creates mol object from SMILES
        if uncharge:
            mol = Chem.MolStandardize.Uncharger().uncharge(mol)  # protonates or deprotonates the mol object
        inchikey = Chem.inchi.MolToInchiKey(mol)  # converts mol object to canonical SMILES
        return inchikey

    @staticmethod
    def remove_isotope_info(_smiles):
        matches = re.findall(r'\[14([CcHh23@]*)\]', _smiles)
        for match in list(set(matches)):
            no_hH = re.sub(r'[hH23@]*', '', match)
            _smiles = _smiles.replace('[14{}]'.format(match), no_hH)
        return _smiles

    @staticmethod
    def remove_stereo_info(_smiles):
        new_smiles = re.sub(r'[\\/]', '', _smiles)
        return new_smiles

    @staticmethod
    def log_transform(input_list):
        log_list = np.log10(input_list)
        log_list[np.isneginf(log_list)] = 0  # replace -inf with 0
        return log_list

    @staticmethod
    def power_10_transform(input_list):
        power_list = np.power(10, input_list)
        return power_list

    # todo: I am temporarily adding some visualizations methods as utils
    # todo: I can move them to visualize after deciding the design

    @staticmethod
    def print_box_plot(df, x, y, tag=''):  # , colors="husl"
        # Error - is_charged
        sns.set(rc={"figure.figsize": (15, 7.5)})
        # sns.plotting_context("paper")
        sns.set_theme(style="whitegrid")
        this = df[[x, y]]
        # sns.boxplot(x=x, y=y, data=this, palette=colors)
        new_y = 'log(DT50)'
        this.rename(columns={y: new_y}, inplace=True)
        sns.boxplot(x=x, y=new_y, data=this)
        plt.savefig(os.path.join(df.export_path, 'output', 'figures', 'boxplot_{}_{}_{}.pdf'.format(x, y, tag)), bbox_inches='tight')
        plt.close()

    @staticmethod
    def warn_disconnection(smiles_list):
        new_smiles_list = []
        for smiles in smiles_list:
            if '.' in str(smiles):
                smiles_split = smiles.split('.')
                smiles_split.sort(key=len, reverse=True)
                cropped_smiles = smiles_split[0]
                print('Warning: There is a disconnection (i.e., ".") in the SMILES notation: {}.'
                      'This can potentially influence calculation of descriptors.'
                      'So only {} will be used'.format(smiles, cropped_smiles))

                new_smiles = cropped_smiles
            else:
                new_smiles = smiles
            new_smiles_list.append(new_smiles)
        return new_smiles_list


    @staticmethod
    def get_display_name(input_str):
        """
        Get a display name ready for output/plotting
        @param input_str: name of descriptor of regressor to be converted to display-ready name
        @return:
        """
        if input_str == 'padel':
            return 'PaDEL'
        elif input_str == 'maccs':
            return 'MACCS'
        elif input_str == 'ep_trig':
            return 'eP_rule'
        elif input_str == 'SGDRegressor':
            return 'SGD'
        elif input_str in ['KNeighborsRegressor', 'KNN Regressor']:
            return 'KNN'
        elif input_str == 'MLPRegressor':
            return 'MLP'
        elif input_str == 'Support Vector Regressor':
            return 'SVR'
        elif input_str in ['GradientBoostingRegressor', 'Gradient Boosting Regressor']:
            return 'GB'
        elif input_str in ['RandomForestRegressor', 'Random Forest Regressor']:
            return 'RF'
        elif input_str == 'AdaBoostRegressor':
            return 'AB'
        elif input_str == 'Gaussian Process Regressor':
            return 'GPR'
        elif input_str == 'pca':
            return 'PCA'
        elif input_str == 'svd':
            return 'SVD'
        elif input_str == 'padel+maccs+ep_trig':
            return 'all'
        else:
            return input_str
