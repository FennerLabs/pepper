from subprocess import Popen, PIPE
import os
import re
import numpy as np
from rdkit import Chem

class Util:
    @staticmethod
    def name_to_smiles(name):
        smiles = name
        return smiles

    @staticmethod
    def canonicalize_smiles(smiles, uncharge=False):  # using rdkit
        mol = Chem.MolFromSmiles(smiles)  # creates mol object from SMILES
        if uncharge:
            mol = Chem.MolStandardize.Uncharger().uncharge(mol)  # protonates or deprotonates the mol object
        new_smiles = Chem.rdmolfiles.MolToSmiles(mol)  # converts mol object to canonical SMILES
        can_smiles = Chem.CanonSmiles(new_smiles)
        return can_smiles

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



