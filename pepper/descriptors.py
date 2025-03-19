import getpass
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import PandasTools

from padelpy import padeldescriptor
padeldescriptor(d_3d=False)

from padelpy import from_smiles
from mordred import Calculator, descriptors

from rdkit.Chem import AllChem as ac
from rdkit.Chem import rdFingerprintGenerator

from sklearn.feature_selection import VarianceThreshold

# sys.path.insert(0, self.get_path_to_enviPath_python() + 'enviPath_python/') # for development only
# sys.path.insert(0, self.get_path_to_enviPath_python())
from enviPath_python import enviPath
from enviPath_python.objects import Package, ParallelCompositeRule, RelativeReasoning

from pepper.pepper import Pepper
from pepper.datastructure import DataStructure


class Descriptors(Pepper):
    def __init__(self, pep: Pepper):
        """
        Initiate Descriptors object
        The descriptors object contains all
        :param pep: Pepper object used to obtain global pepper settings
        """
        super().__init__()
        self.pepper = pep
        self.set_data_directory(os.path.join(self.pepper.data_directory, 'descriptors'))
        self.model_data = pd.DataFrame()

        # attributes from pepper
        self.data_type = pep.data_type
        self.tag = pep.tag
        self.id_name = pep.id_name
        self.smiles_name = pep.smiles_name

        # Storing data for the different types of descriptors
        self.maccs = pd.DataFrame()
        self.maccs_tsv = self.build_output_filename('MACCS')
        self.padel = pd.DataFrame()
        self.padel_tsv = self.build_output_filename('PaDEL')
        self.mordred = pd.DataFrame()
        self.mordred_tsv = self.build_output_filename('mordred')
        self.ep_trig = pd.DataFrame()
        self.ep_trig_tsv = self.build_output_filename('enviPath_triggered')
        self.ep_prob = pd.DataFrame()
        self.ep_prob_tsv = self.build_output_filename('enviPath_probability')
        self.koc = pd.DataFrame()
        self.koc_tsv = self.build_output_filename('Koc')
        self.mfps = pd.DataFrame()
        self.mfps_tsv = self.build_output_filename('mfps')
        self.clogp = pd.DataFrame()
        self.clogp_tsv = self.build_output_filename('clogp')
        self.plant_fp = pd.DataFrame()
        self.plant_fp_tsv = self.build_output_filename('plant_fp')
        self.rdkitfps = pd.DataFrame()
        self.rdkitfps_tsv = self.build_output_filename('RDKit_fps')

        self.features = pd.DataFrame()


        self.feature_space_list = ['maccs', 'padel', 'ep_trig',
                              'ep_prob', 'mfps', 'mordred',
                              'clogp', 'plant_fp', 'rdkitfps']
        self.current_feature_space = '' # e.g., 'maccs', 'padel+ep_trig', 'all'
        self.feature_space_map = {} # links feature names to feature space for visualisation e.g., {'struct-1: 'maccs'}

        # enviPath settings
        #####################
        # Default package for triggered rules: EAWAG-BBD
        self.ep_trig_rule_package = 'http://envipath.org/package/32de3cf4-e3e6-4168-956e-32fa5ddb0ce1'
        # Alternatively, the following package can be used. It includes the newer soil-specific rules.
        # For this, access and login is currently required:
        # self.ep_trig_rule_package = 'https://envipath.org/package/55fa3a97-db19-442f-8108-954f7be95e1c'

        # Default model for rule probability calculation: BBD - ECC - Multi - 2023-09-05
        self.ep_prob_relative_reasoning_id = 'https://envipath.org/package/32de3cf4-e3e6-4168-956e-32fa5ddb0ce1/' \
                                             'relative-reasoning/23e1b2ec-dcc0-4389-9b65-afd52bd72e27'

        # Default enviPath instance
        self.instance_host = 'https://envipath.org/'

        # For the default settings, no login to enviPath is required. However, if the package or relative reasoning
        # for ep_trig or ep_prob is not public, then this ep_login_required should be set to True
        self.ep_login_required = False

    def set_data(self, data: DataStructure):
        """
        Set DataStructure object for which descriptors should be calculated
        :param data: DataStructure object used to obtain model data frame and data type
        """
        self.model_data = data.get_model_data()
        self.data_type = data.get_data_type()
        # self.name_data = data.get_name_data()
        self.id_name = data.id_name
        self.smiles_name = data.smiles_name


    def load_descriptors(self, from_csv=False,
                         PaDEL=False, mordred=False,
                         MACCS=False, enviPath_prob=False, enviPath_trig=False,
                         mfps=False, RDKit_fps=False,
                         plant_fingerprints=False,
                         Koc=False,  clogp=False,
                         load_by_feature_name = False,
                         feature_name_list: list = None,
                         feature_space_map: dict = None):
        """
        Load the desired descriptors
        :param from_csv: If False, calculates them instead of loading from csv. The first run must be False.
        :param PaDEL: if True, load from csv file or calculate
        :param MACCS: if True, load from csv file or calculate
        :param enviPath_prob: if True, load from csv file or calculate
        :param enviPath_trig: if True, load from csv file or calculate
        :param mordred: if True, load from csv file or calculate
        :param Koc: if True, load from provided csv file
        :param mfps: if True, apply Morgan (m) Fingerprints (fps) with radius = 2
        :param clogp: if True, retrieve cLogP from database (originally calculated using OPERA)
        :param RDKit_fps: if True, calculate rdkit fingerprints
        :param plant_fingerprints: if True, include fingerprints to describe plants
        :param load_by_feature_name: list of feature names to be loaded. If provided, only specified features are calculated
        :param feature_space_map: feature space map from another Descriptors object
        """
        print("\n############# Load descriptors ############# ")

        # if specific features are provided, they are mapped to feature spaces (e.g., 'maccs', 'padel')
        if load_by_feature_name:
            print("Calculate features by feature name")
            features_to_be_calculated = {}  #sort features by feature_space to provide to load functions
            for feature_name in feature_name_list:
                if feature_space_map[feature_name] not in features_to_be_calculated.keys():
                    features_to_be_calculated[feature_space_map[feature_name]] = []
                features_to_be_calculated[feature_space_map[feature_name]].append(feature_name)
            [print('{}: {}'.format(f, len(features_to_be_calculated[f]))) for f in features_to_be_calculated.keys()]
            # set features spaces to be covered
            if 'maccs' in features_to_be_calculated.keys():
                MACCS = True
            if 'padel' in features_to_be_calculated.keys():
                PaDEL = True
            if 'ep_trig' in features_to_be_calculated.keys():
                enviPath_trig = True
            if 'ep_prob' in features_to_be_calculated.keys():
                enviPath_prob = True
            if 'mfps' in features_to_be_calculated.keys():
                mfps = True
            if 'clogp' in features_to_be_calculated.keys():
                clogp = True
            if 'plant_fp' in features_to_be_calculated.keys():
                plant_fingerprints = True
            if 'Koc' in features_to_be_calculated.keys():
                Koc = True
            if 'mordred' in features_to_be_calculated.keys():
                mordred = True
            if 'rdkitfps' in features_to_be_calculated.keys():
                RDKit_fps = True

        if MACCS: # maccs fingerprints are always calculated as a whole matrix
            self.maccs_tsv = self.build_output_filename('MACCS')
            if from_csv:
                self.maccs = pd.read_csv(self.maccs_tsv, sep='\t')
            elif not self.model_data.empty:
                self.calculate_MACCS_fingerprints()
            else:
                raise ValueError('No model data nor csv file provided')
            self.populate_feature_space_map('maccs', self.maccs.columns.values)

        if PaDEL:
            if from_csv:
                self.padel = pd.read_csv(self.padel_tsv, sep='\t')
                # temporary workaround to allow loading precalculated descriptor files without the "PaDEL" tag
                # todo: remove workaround
                if 'PaDEL' not in self.padel.columns[0]:
                    self.padel.columns = [f'PaDEL-{column_name}' for column_name in self.padel.columns]
                    self.padel.rename(columns = {'PaDEL-' + self.smiles_name : self.smiles_name}, inplace=True)
            elif not self.model_data.empty:
                self.calculate_PaDEL_descriptors()
            else:
                raise ValueError('No model data nor csv file provided')
            self.populate_feature_space_map('padel', self.padel.columns.values)

        if mordred:
            self.mordred_tsv = self.build_output_filename('mordred')
            if from_csv:
                self.mordred = pd.read_csv(self.mordred_tsv, sep='\t')
            elif load_by_feature_name and not self.model_data.empty:
                self.calculate_mordred_descriptors(features_to_be_calculated['mordred'])
            elif not self.model_data.empty:
                self.calculate_mordred_descriptors()
            else:
                raise ValueError('No model data nor csv file provided')
            self.populate_feature_space_map('mordred', self.mordred.columns.values)

        if enviPath_trig:
            if from_csv:
                self.ep_trig_tsv = self.build_output_filename('enviPath_triggered')
                self.ep_trig = pd.read_csv(self.ep_trig_tsv, sep='\t', index_col=False)
            elif load_by_feature_name and not self.model_data.empty:
                self.calculate_enviPath_descriptors(triggered=True,
                                                    feature_name_list = features_to_be_calculated['ep_trig'])
            elif not self.model_data.empty:
                self.calculate_enviPath_descriptors(triggered=True)
            else:
                raise ValueError('No model data nor csv file provided')
            self.populate_feature_space_map('ep_trig', self.ep_trig.columns.values)

        if enviPath_prob:
            if from_csv:
                self.ep_prob_tsv = self.build_output_filename('enviPath_probability')
                self.ep_prob = pd.read_csv(self.ep_prob_tsv, sep='\t')
            elif load_by_feature_name and not self.model_data.empty:
                self.calculate_enviPath_descriptors(probabilities=True,
                                                    feature_name_list = features_to_be_calculated['ep_prob'])
            elif not self.model_data.empty:
                self.calculate_enviPath_descriptors(probabilities=True)
            else:
                raise ValueError('No model data nor csv file provided')
            self.populate_feature_space_map('ep_prob', self.ep_prob.columns.values)

        if Koc:
            if from_csv:
                self.koc_tsv = self.build_output_filename('Koc')
                self.koc = pd.read_csv(self.koc_tsv, sep='\t')  # todo: add method to "clean" the opera output
            else:
                raise NotImplementedError('There is currently no implementation to estimate Koc automatically'
                                          'Please provide Koc values (from, e.g., OPERA) in a standard descriptor '
                                          'file format')
            self.populate_feature_space_map('koc', self.koc.columns.values)

        if mfps:
            if from_csv:
                self.mfps_tsv = self.build_output_filename('mfps')
                self.mfps = pd.read_csv(self.mfps_tsv, sep='\t')
            elif load_by_feature_name and not self.model_data.empty:
                self.calculate_morgan_fingerprints(features_to_be_calculated['mfps'])
            elif not self.model_data.empty:
                self.calculate_morgan_fingerprints()
            else:
                raise ValueError('No model data nor csv file provided')
            self.populate_feature_space_map('mfps', self.mfps.columns.values)

        if RDKit_fps:
            if from_csv:
                self.rdkitfps_tsv = self.build_output_filename('RDKit_fps')
                self.rdkitfps = pd.read_csv(self.rdkitfps_tsv, sep='\t')
            elif load_by_feature_name and not self.model_data.empty:
                self.calculate_rdkit_fingerprints(features_to_be_calculated['RDKit_fps'])
            elif not self.model_data.empty:
                self.calculate_rdkit_fingerprints()
            else:
                raise ValueError('No model data nor csv file provided')
            self.populate_feature_space_map('rdkitfps', self.rdkitfps.columns.values)

        if clogp:
            if from_csv:
                self.clogp_tsv = self.build_output_filename('clogp')
                self.clogp = pd.read_csv(self.clogp_tsv, sep='\t')
            elif not self.model_data.empty:
                self.calculate_clogp()
            else:
                raise ValueError('No model data nor csv file provided')
            self.populate_feature_space_map('clogp', self.clogp.columns.values)

        if plant_fingerprints:
            try:
                plant_fp_df = pd.read_csv(self.plant_fp_tsv, sep='\t')
                plant_fp_df.drop(columns=['dataset'], inplace=True)
                self.plant_fp = pd.merge(self.model_data['plant'],
                                         # plant_fp_df[['pca0', 'pca1', 'plant']],
                                         plant_fp_df,
                                         how='left', on='plant')
                print(self.plant_fp.shape)
            except ValueError:
                print('Could not find plant_fp file')
            self.populate_feature_space_map('plant_fp', self.plant_fp.columns.values)

        return

    def merge_descriptors(self, descriptors_list: list = None):
        """
        Merge descriptors into features matrix to be used for modelling
        :param descriptors_list: List of descriptor types to be merged.
        If None then the method creates a list based on calculated descriptors;
        descriptors_list is None when feature space is 'all'.
        """
        if descriptors_list is None:
            descriptors_list = []
            if not self.maccs.empty:
                descriptors_list.append('maccs')
            if not self.padel.empty:
                descriptors_list.append('padel')
            if not self.mordred.empty:
                descriptors_list.append('mordred')
            if not self.ep_trig.empty:
                descriptors_list.append('ep_trig')
            if not self.ep_prob.empty:
                descriptors_list.append('ep_prob')
            if not self.koc.empty:
                descriptors_list.append('koc')
            if not self.mfps.empty:
                descriptors_list.append('mfps')
            if not self.clogp.empty:
                descriptors_list.append('clogp')
            if not self.plant_fp.empty:
                descriptors_list.append('plant_fp')
            if not self.rdkitfps.empty:
                descriptors_list.append('rdkitfps')

        self.features = self.get_features_by_keyword(descriptors_list[0])
        if len(descriptors_list) > 1:
            for k in descriptors_list[1:]:
                self.features = self.features.merge(self.get_features_by_keyword(k), on=self.smiles_name)


    def get_features_by_keyword(self, keyword):
        """
        Get the descriptor table by descriptor type (keyword)
        :param keyword: descriptor type
        :return: Descriptors of the indicated descriptor type (keyword) as a DataFrame
        """

        if keyword == 'maccs':
            feature_matrix = self.maccs
        elif keyword == 'padel':
            feature_matrix = self.padel
        elif keyword == 'mordred':
            feature_matrix = self.mordred
        elif keyword == 'ep_trig':
            feature_matrix = self.ep_trig
        elif keyword == 'ep_prob':
            feature_matrix = self.ep_prob
        elif keyword == 'plant_fp':
            feature_matrix = self.plant_fp
        elif keyword == 'koc':
            return self.koc
        elif keyword == 'mfps':
            return self.mfps
        elif keyword == 'clogp':
            return self.clogp
        elif keyword == 'rdkitfps':
            return self.rdkitfps
        else:
            raise NotImplementedError('No features available for {}'.format(keyword))
        assert not feature_matrix.empty, ('The feature space {} is currently not loaded'.format(keyword))
        return feature_matrix

    def check_if_loaded(self, feature: str):
        """
        Check if a given feature is loaded to the Descriptors object
        :param feature: string of feature to be che checked, e.g., "maccs"
        :return: True if loaded, False otherwise
        """
        if self.get_features_by_keyword(feature).empty:
            return False
        else:
            return True

    def define_feature_space(self, keyword: str):
        """
        Define type of feature space to create from available descriptors
        :param keyword: type of feature space - 'all', 'maccs', 'padel', 'mordred' ,'ep_trig', 'ep_prob', 'mfps', clogp', or any combination of them (e.g., 'maccs+mordred")
        """
        print("Define feature space:", keyword)
        self.current_feature_space = keyword

        if keyword == 'all':
            self.merge_descriptors()
        elif keyword in self.feature_space_list:
            self.features = self.get_features_by_keyword(keyword)
        elif '+' in keyword:
            split = keyword.split('+')
            split_checked = []
            for desc in split:
                if desc in self.feature_space_list:
                    split_checked.append(desc)
                else:
                    raise Warning('Feature space keyword "{}" not found in possible feature space'.format(desc))
            self.merge_descriptors(split_checked)
        else:
            raise ValueError("The provided feature space keyword does not exist. It should be one of the following, "
                             "or a combination of concatenated strings with '+' as a separator:"
                             " 'all', 'maccs', 'padel', 'mordred' ,'ep_trig', 'ep_prob', 'clogp', 'mfps', 'rdkitfps")

        # # Define the model feature space by taking a subset of descriptors.features based on model data.
        return

    def get_current_feature_space(self):
        """
        Get currently used feature space
        :return: str describing feature space
        """
        return self.current_feature_space

    def get_feature_space_list(self):
        return self.feature_space_list

    def calculate_MACCS_fingerprints(self):
        """
        Calculate MACCS fingerprints via RDKit
        For interpretation (SMARTS patterns corresponding to fingerprints), see here:
        https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/MACCSkeys.py
        """
        print('-> calculate MACCS fingerprints')
        # format header
        header = []
        fp = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(np.array(self.model_data[self.smiles_name])[0])).ToList()
        [header.append('struct-{}'.format(i + 1)) for i, e in enumerate(fp[1:])]
        D = {}
        for index, row in self.model_data.iterrows():
            smiles = row[self.smiles_name]

            try:
                fp = Chem.MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
                descriptors_id = fp.ToList()[1:]
            except RuntimeError:
                print('Warning: No MACCS descriptor could be calculated for compound {}, smiles = {}'.format(smiles, row[
                    'reduced_smiles']))
            else:
                desc_D = {}
                for i, h in enumerate(header):
                    desc_D[h] = descriptors_id[i]
                D[smiles] = desc_D  # here the keys of the dictionary are set.

        self.maccs = pd.DataFrame.from_dict(D, columns=header, orient='index')
        self.maccs[self.smiles_name] = self.maccs.index
        self.maccs.to_csv(self.maccs_tsv, sep='\t', index=False)

        return

    def calculate_PaDEL_descriptors(self):
        """
        Calculate PaDEL descriptors via padelpy (DOI: 10.1002/jcc.21707)
        """
        print('-> calculate PaDEL descriptors')
        D = {}
        smiles_list = []
        for index, row in self.model_data.iterrows():
            smiles = row[self.smiles_name]
            try:
                padel_descriptors = from_smiles(smiles)
            except RuntimeError:
                print('Warning: No PaDEL descriptor could be '
                      'calculated for compound: {}'.format(smiles))
            else:
                D[smiles] = padel_descriptors
                smiles_list.append(smiles)

        self.padel = pd.DataFrame.from_dict(D, orient='index')
        self.padel.columns = [f'PaDEL-{column_name}' for column_name in self.padel.columns]
        self.padel[self.smiles_name] = smiles_list
        self.padel.to_csv(self.padel_tsv, sep='\t', index=False)

    def calculate_enviPath_descriptors(self, triggered=False, probabilities=False, feature_name_list = None):
        """
        Obtain descriptors from enviPath via enviPath-python
        :param triggered: if True, ep_trig (triggered rules) is calculated
        :param probabilities: if True, ep_prob (rule probabilities) is calculated
        :param feature_name_list
        """
        print('-> calculate enviPath rule descriptors')
        eP = enviPath(self.instance_host)

        # logging in to envipath
        if self.ep_login_required:
            user = input('Username for envipath.org: ')
            password = getpass.getpass()
            eP.login(user, password)

        # Load packages BBD and new soil rules
        if triggered:
            print("Calculating enviPath triggered rules...")
            eP_package = Package(eP.requester, id=self.ep_trig_rule_package)
            rule_ids = self.get_composite_rules(eP_package)
            rules_list = []
            for r_id in rule_ids:
                rule = ParallelCompositeRule(eP.requester, id=r_id)
                name = rule.get_name() + '-trig'
                if feature_name_list is not None and name not in feature_name_list:
                    continue #skip features that are not in list
                rules_list.append(rule)
            self.get_rule_descriptors(rules_list)

        # load relative reasoning for prob values
        if probabilities:
            print("Calculating enviPath rule probabilities...")
            relres = RelativeReasoning(eP.requester, id=self.ep_prob_relative_reasoning_id)
            self.get_rule_probabilities(relres, feature_name_list)

    @staticmethod
    def get_composite_rules(package, list_of_rules=None):
        """
        Extract the bt rules from the enviPath data package
        :param package:
        :param list_of_rules:
        :return:
        """
        rule_ids = []
        rules = package.get_rules()
        for r in rules:
            if list_of_rules:
                if r.get_name() in list_of_rules:
                    rule_ids.append(r.get_id())
            else:
                if r.is_composite_rule():
                    rule_ids.append(r.get_id())
        return rule_ids

    def get_rule_descriptors(self, list_of_rules):
        """
        Applies each bt rules to all the SMILES in model data, and obtain a boolean matrix of triggered rules
        :param list_of_rules: list of enviPath ParallelRule objects
        """
        D = {}  # D{'rule': {'ID' : value}, ...}

        for rule in list_of_rules:
            print(rule)
            name = rule.get_name()
            D[name + '-trig'] = {}

            smiles_list = []
            for index, row in self.model_data.iterrows():
                smiles = row[self.smiles_name]
                try:
                    out = rule.apply_to_smiles(smiles)
                except:
                    print('Could not process SMILES:', smiles)
                    #break
                if out == []:
                    value = 0
                else:
                    value = 1
                D[name+'-trig'][index] = value
                smiles_list.append(smiles)
        self.ep_trig = pd.DataFrame(D)
        self.ep_trig[self.smiles_name] = smiles_list
        self.ep_trig.to_csv(self.ep_trig_tsv, sep='\t', index=False)

    def get_rule_probabilities(self, relative_reasoning, feature_list_prob):
        """
        #todo: this is porably not working
        Expand each SMILES in the model_data using the relative_reasoning, and get the reaction probability
        :param relative_reasoning: relative reasoning model to be used to obtain probabilities
        """
        D = {}
        rules = set([])
        smiles_list = []
        for index, row in self.model_data.iterrows():
            smiles = row[self.smiles_name]
            # try:
            out = self.expand_smiles(smiles, relative_reasoning)
            # except:
            #     print('Could not process SMILES:', smiles)
            # else:
            res_dict = self.result_to_rules_dict(out, feature_list_prob)
            for rule in res_dict.keys():
                rules.add(rule)
            D[index] = res_dict
            smiles_list.append(smiles)
        self.ep_prob = pd.DataFrame(columns=list(rules).sort())  # todo: check this sort

        for ID in D.keys():
            for rule_name in D[ID].keys():
                probability = D[ID][rule_name]
                row_index = self.ep_prob[self.ep_prob[self.id_name] == ID].index
                self.ep_prob.loc[row_index, rule_name] = probability

        self.ep_prob[self.smiles_name] = smiles_list
        self.ep_prob.fillna(0, inplace=True)
        self.ep_prob.to_csv(self.ep_prob_tsv, sep='\t', index=False)

    def calculate_mordred_descriptors(self):
        """"
        Calculate all molecular descriptors available in mordred
        Parameters:
        dataframe: A data frame containing SMILES
        ------------------------------------------------------
        in mordred, pandas means calculate descriptors from a mol object
        """
        print('-> calculate mordred descriptors')
        dataframe = self.model_data
        if 'ROMol' not in dataframe.columns:
            PandasTools.AddMoleculeColumnToFrame(dataframe, self.smiles_name, 'ROMol', includeFingerprints=True)
            print("calculated ROMol")
        else:
            pass

        calc_all = Calculator(descriptors, ignore_3D=True)
        mordred_descriptors = pd.DataFrame(calc_all.pandas(dataframe['ROMol']))
        self.mordred = mordred_descriptors
        self.mordred.columns = [f'mordred-{column_name}' for column_name in self.mordred.columns]
        # Include CanonicalSMILES
        self.mordred[self.smiles_name] = self.model_data[self.smiles_name]
        self.mordred.to_csv(self.mordred_tsv, sep='\t', index=False)

    def calculate_morgan_fingerprints(self):
        """
        Calculates the Morgan fingerprints (as ECFPs, because radius = 2) from RDKit
        for each compound.
        :return: a dataframe file with the morgan fingerprints of all compounds.
        """
        print('->calculate morgan fingerprints')
        dataframe = self.model_data
        if 'mfps' not in dataframe.columns:
            smiles_list = dataframe[self.smiles_name]
            # Converting standardised SMILES to Morgan fingerprints (ECFP4s)
            mf_bv = []  # morgan fingerprints as BitVectors
            for smiles in smiles_list:
                mol = Chem.MolFromSmiles(smiles)
                mf_bv.append(ac.GetMorganFingerprintAsBitVect(mol, 2, 2048))
            mfps = np.array([list(x) for x in mf_bv])
            print("calculated morgan fingerprints")
        else:
            pass
        mfps_descriptors = pd.DataFrame(mfps)
        self.mfps.columns = [f'Morgan-{column_name}' for column_name in self.mfps.columns]  # todo: Add names to columns otherwise they collide
        mfps_descriptors[self.smiles_name] = dataframe[self.smiles_name]
        self.mfps = mfps_descriptors
        self.mfps.to_csv(self.mfps_tsv, sep='\t', index=False)

    def calculate_rdkit_fingerprints(self):
        """
        Calculates the Morgan fingerprints (as ECFPs, because radius = 2) from RDKit
        for each compound.
        :return: a dataframe file with the morgan fingerprints of all compounds.
        """
        print('->calculate rdkit fingerprints')
        rdkgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)

        dataframe = self.model_data
        smiles_list = dataframe[self.smiles_name]

        # Converting standardised SMILES to Morgan fingerprints (ECFP4s)
        rdkit_CountVector_list = []  # rdkit fingerprints as CountVectors
        rdkit_BitVector_list = []  # morgan fingerprints as BitVectors
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            # rdkit_CountVector_list.append(rdkgen.GetCountFingerprintAsNumPy(mol))
            rdkit_BitVector_list.append(rdkgen.GetFingerprintAsNumPy(mol))
        print("calculated rdkit fingerprints")

        # self.rdkitfps = pd.DataFrame(rdkit_CountVector_list)
        self.rdkitfps = pd.DataFrame(rdkit_BitVector_list)
        self.rdkitfps.columns = [f'RDKit-{column_name}' for column_name in self.rdkitfps.columns]
        self.rdkitfps[self.smiles_name] = dataframe[self.smiles_name]
        self.rdkitfps.to_csv(self.rdkitfps_tsv, sep='\t', index=False)
        return

    def calculate_clogp(self):
        """
        Retrieves the CLogP from database for each compound.
        :return: a dataframe file with the CLogP descriptors of all compounds.
        """
        print('->get CLogP') #todo: remove path to user data (@jose)
        clogp_database = pd.read_csv('/Users/corderjo/pepper_data/data_structure'
                                     '/name_data_ActivatedSludge_WWTP_CompleteDatabase.tsv',
                                     sep='\t')[[self.smiles_name, 'ClogP']]

        self.clogp = clogp_database
        self.clogp.to_csv(self.clogp_tsv, sep='\t', index=False)

    @staticmethod
    def expand_smiles(smiles, rr):
        """
        Get all potential TPs by applying enviPath biotransformation and relative reasoning rules
        :param smiles: input smiles
        :param rr: relative reasoning object
        :return: list of dictionaries for each predicted TP: {'smiles': smiles,
                                                              'name': rule name,
                                                              'probability': relative reasoning probability}
        """
        res = rr.classify_smiles(smiles)
        # sort by probability
        res.sort(reverse=True, key=lambda x: x['probability'])
        return res

    @staticmethod
    def result_to_rules_dict(result, feature_list_prob):
        """
        Translates result from enviPath node expansion into a compound dictionary
        :param result: list of dictionaries with predicted TP information
        :return: dictionary of TPs
        """
        rules_dict = {}
        for r in result:
            probability = float(r['probability'])
            name = r['name'] + '-prob'
            if feature_list_prob is not None and name not in feature_list_prob:
                continue  # skip features that are not in list
            if name not in rules_dict.keys():
                rules_dict[name] = probability
            else:
                if probability > rules_dict[name]:
                    rules_dict[name] = probability
        return rules_dict

    def simple_preprocessing(self):
        print("-> Removing features with zero variance (still in preparation)")
        print("Shape before removing zero-variance:{}".format(self.features.shape))
        # VarianceThreshold().fit_transform(np.array(self.features))
        # print("Shape after removing zero-variance:{}".format(self.features.shape))

    def populate_feature_space_map(self, feature_space, list_of_feature_names):
        for feature_name in list_of_feature_names:
            if feature_name in [self.id_name, self.smiles_name]: # do not consider smiles and id column headers
                continue
            assert feature_name not in self.feature_space_map, f"Duplicate feature name in {feature_space} feature space, please check {feature_name}"
            self.feature_space_map[feature_name] = feature_space

        return