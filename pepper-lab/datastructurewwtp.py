import os
import copy
import matplotlib.pyplot as plt
import yaml

from pepper.pepper import Pepper
from pepper.util import Util
from pepper.datastructure import DataStructure


import pandas as pd  # Do we need to keep importing these general packages?
import numpy as np
import pubchempy as pcp


# Some variables
plants_group_1 = ['S1', 'S9', 'S29', 'S51']
plants_group_2 = ['S2', 'S5', 'S11', 'S24', 'S25', 'S40']
plants_group_3 = ['S10', 'S39', 'S53']
long_HRT_plant_list = ['S66', 'S67']
plant_list = plants_group_1 + plants_group_2 + plants_group_3 + long_HRT_plant_list
amar_AS_plant_list = ['lug', 'wer', 'neu', 'enk', 'kap', 'kat']


class DataStructureWWTP(DataStructure):
    def __init__(self, pep: Pepper):
        super().__init__(pep)
        self.data_type = 'WWTP'
        self.curation_type = pep.get_curation_type()
        self.setup_name = pep.setup_name
        self.compound_name = pep.get_compound_name()
        self.id_name = pep.get_id_name()
        self.smiles_name = pep.get_smiles_name()
        self.plant_name = 'plant'
        self.retention_time_name = 'RT [min]'
        self.LOQ_flag_name = 'LOQasEffluent'

        # The difference with respect to DataStructure is to allow having different names according to the config.
        # This can mean different set of plants and thus different data except for the raw data.
        self.raw_data_tsv = self.build_output_filename('raw_data')
        self.full_data_tsv = self.build_output_filename('full_data_' + self.setup_name)
        self.cpd_data_tsv = self.build_output_filename('cpd_data_' + self.setup_name)
        self.model_data_tsv = self.build_output_filename('model_data_' + self.setup_name + '_' + self.curation_type)
        self.name_data_tsv = self.build_output_filename('name_data_' + self.setup_name)
        self.joint_data_tsv = self.build_output_filename('joint_data_' + self.setup_name)
        self.opera_properties_tsv = self.build_output_filename('opera_properties')

        # WWTP specific attributes
        self.name_data = pd.DataFrame()
        self.joint_data = pd.DataFrame()
        self.opera_properties = pd.DataFrame()

        self.single_model_data = pd.DataFrame()
        self.single_model_data_tsv = self.build_output_filename('single_model_data')

        self.one_plant_df = pd.DataFrame()

        self.duplicates_df = pd.DataFrame()

        # choosing plants to work with
        self.plant_list = []
        self.ndn_plants = ['S2', 'S5', 'S11', 'S24', 'S40',
                           'Katrineholm', 'Käppala', 'Luggage Point', 'Neugut', 'Werdhölzli',
                           'alt', 'bir', 'ehr', 'kol',
                           'WWTP_1', 'WWTP_2', 'WWTP_3', 'WWTP_4', 'WWTP_5',
                           'WWTP_6', 'WWTP_7', 'WWTP_8', 'WWTP_9', 'WWTP_10',
                           'WWTP_11', 'WWTP_12', 'WWTP_13', 'WWTP_14', 'WWTP_15',
                           'WWTP_16', 'WWTP_17', 'WWTP_18', 'WWTP_19', 'WWTP_20',
                           'WWTP_21', 'WWTP_22', 'WWTP_23', 'WWTP_24', 'WWTP_25',
                           'WWTP_26']

        self.c_eliminating_plants = ['S25', 'Enköping', 'air', 'oli']

    def set_plant_list(self, my_plant_list=False):
        """
        Only applicable to wastewater treatment plant related data.
        Allows to select a subset of the whole database based on a list of treatment plants of interest.
        The defaults behaviour creates a list of plants that are considered to be similar
        in terms of technology and performance.
        :param my_plant_list:
        """
        if my_plant_list:
            self.plant_list = my_plant_list
        else:
            if self.tag == 'aus_data':
                self.plant_list = ['S1', 'S10', 'S11', 'S2', 'S24', 'S25',
                                   'S29', 'S39', 'S40', 'S5', 'S51', 'S53', 'S9']
            elif self.tag == 'amar_data':
                self.plant_list = ['Enköping', 'Käppala', 'Katrineholm', 'Werdhölzli', 'Neugut', 'Luggage Point']
            elif self.tag == 'galpro_data':
                self.plant_list = ['Acino', 'Amino', 'Akorn', 'Sintetica', 'SwissCaps', 'SwissCo', 'Corden']
            elif self.tag == 'snf_data':
                self.plant_list = ['alt', 'bir', 'ehr', 'kol', 'oli', 'air']
            elif self.tag == 'swe2_data':
                self.plant_list = list(self.raw_data[self.plant_name].unique())
            else:
                print("unknown default plant list OR unrecognized tag. "
                      "Allowed tags are 'aus_data', 'amar_data', 'snf_data' and 'swe2_data'.")

    def get_plant_list(self):
        self.set_plant_list()
        return self.plant_list

    def get_target_variable_from_inf(self, target_variable_list=None, verbose=False):

        # Avoid errors from missing data
        # Replace n.d. with NaN in columns A and B
        self.full_data['influent'].replace('n.d.', np.nan, inplace=True)
        self.full_data['effluent'].replace('n.d.', np.nan, inplace=True)

        # Convert columns A and B to numeric types
        self.full_data['influent'] = pd.to_numeric(self.full_data['influent'], errors='coerce')
        self.full_data['effluent'] = pd.to_numeric(self.full_data['effluent'], errors='coerce')

        if verbose:
            print("Calculating endpoints from influent and effluent values")

        if target_variable_list is None:
            target_variable_list = self.target_variable_list

        if 'B' in target_variable_list:
            self.full_data['B'] = self.full_data['effluent'] / self.full_data['influent']

        if 'B(%)' in target_variable_list:
            self.full_data['B(%)'] = (self.full_data['effluent'] / self.full_data['influent']) * 100

        if 'R(%)' in target_variable_list:
            self.full_data['R(%)'] = ((self.full_data['influent'] - self.full_data['effluent'])
                                      / self.full_data['influent']) * 100

        if 'logB(%)' in target_variable_list:  # notice this is the log base 10 of the percentage
            self.full_data['logB(%)'] = np.log10((self.full_data['effluent'] / self.full_data['influent']) * 100)

        if 'logB' in target_variable_list:  # notice this is the log base 10 of the percentage
            self.full_data['logB'] = np.log10((self.full_data['effluent'] / self.full_data['influent']))

        valid_endpoints = {'B', 'B(%)', 'R(%)', 'logB', 'logB(%)'}
        if not set(target_variable_list).intersection(valid_endpoints):
            raise ValueError(f"No endpoint was calculated: valid endpoints are {', '.join(valid_endpoints)}")

        return

    def name_to_smiles(self, ignore_not_found=False, verbose=False):
        """Adds (or updates) a SMILES column in name_data based on the SMILES entry in PubChem.
        The names of compounds for which a PubChem entry was not found are displayed too
        """
        if verbose:
            print("\n############# Getting SMILES from PubChem ############# ")

        name_list = self.name_data[self.compound_name]

        smiles_list = []
        inchikey_list = []

        not_found = []

        for name in name_list:
            pcp_properties = pcp.get_properties(['IsomericSMILES', 'InChIKey'], name, 'name', as_dataframe=False)
            if not pcp_properties:
                pcp_properties = pcp.get_properties(['IsomericSMILES', 'InChIKey'], name,
                                                    'inchikey', as_dataframe=False)
                if not pcp_properties:
                    not_found.append(name)

                    pcp_properties = [{'IsomericSMILES': np.nan,
                                       'InChIKey': np.nan}]

            smiles = pcp_properties[0].get('IsomericSMILES')
            inchikey = pcp_properties[0].get('InChIKey')

            smiles_list.append(smiles)
            inchikey_list.append(inchikey)

        self.name_data[self.smiles_name] = smiles_list
        self.name_data[self.inchikey_name] = inchikey_list

        if verbose:
            print("names not found in PubChem:")
            print(not_found)

        if ignore_not_found:
            self.name_data.dropna(subset=self.inchikey_name, inplace=True)
            return

    # These are AUS specific functions

    def curate_annotate(self, load_chemstruct: bool = False, source: str = None, verbose=False):
        """
        New method to curate raw data.
        It is divided in two major methods:
        1) curate by chemical structure, which intends to handle problems that may arise from disconnections in SMILEs,
        wrong annotations, presence of duplicates after simplifying chemicals in terms of stereochemistry.
        2) curate by target variable values, which intends to, flag compounds suspected to form during
        wastewater treatment, deal with extreme values, flag compounds with two datapoints or less.
        :param verbose: If True provides detailed messages of the curation steps
        :param source: If csv loads from csv....
        :param load_chemstruct: If True, curated smiles are loaded from csv. The first run it must be False
        :return: None
        """
        if source == 'pepper_data':
            assert os.path.exists(self.full_data_tsv), "Error: file {} does not exist".format(self.full_data_tsv)
            self.full_data = pd.read_csv(self.full_data_tsv, sep='\t', )
            assert os.path.exists(self.name_data_tsv), "Error: file {} does not exist".format(self.name_data_tsv)
            self.name_data = pd.read_csv(self.name_data_tsv, sep='\t')
            print('Existing curated files loaded from {} and {}'.format(self.full_data_tsv, self.name_data_tsv))
            return

        elif not source:

            self.get_plant_list()
            self.pre_curation()

            self.curate_by_chemstruct(from_csv=load_chemstruct, verbose=verbose)
            self.set_id_name('Combined_ID')
            self.set_smiles_name('CanonicalSMILES')

            # Get the expected target variables for WWTP data
            if 'influent' in self.full_data.columns:
                self.get_target_variable_from_inf(['B'])
            else:
                if self.tag != 'swe2_data':
                    raise ValueError('influent not found in data')
                else:
                    print("In the case of swe2_data (Yijing's data) we just trust the calculated Breakthrough")
                    pass

            self.merge_duplicates()
            self.curate_by_target_variable(verbose=verbose)

            print("\n############# Saving outputs ############# ")

            if verbose:
                print("\n############# Get Calculated logP from RDKit Crippen module ############# ")
            self.name_data['ClogP'] = self.get_rdkit_properties(clogp=True)

            # save
            self.full_data.to_csv(self.full_data_tsv, sep='\t', index=False)
            self.name_data.to_csv(self.name_data_tsv, sep='\t', index=False)
            self.cpd_data.to_csv(self.cpd_data_tsv, sep='\t', index=False)

        else:
            raise ValueError('source unknown')

        return

    def pre_curation(self):
        """ The purpose of this method is to document how each one goes from the raw data
        ,as shared by collaborators, to a shape compatible with PEPPER's workflow.   """
        self.name_data = self.full_data

        if self.tag == 'aus_data':
            self.shape_aus_full_data()
            self.full_data[self.LOQ_flag_name] = False

        elif self.tag == 'amar_data':
            self.format_amar_data()
            self.correct_by_LODQ(verbose=True)

        elif self.tag == 'snf_data':
            self.format_snf_data()
            self.snf_correct_by_LODQ(verbose=True)

        elif self.tag == 'swe2_data':
            self.compound_name = 'Compound'
            self.target_variable_name = 'B'
            self.get_info_from_raw_data()

        else:
            self.get_info_from_raw_data()

        self.full_data['dataset'] = self.tag
        return

    def get_info_from_raw_data(self):
        try:
            self.name_data = self.raw_data.loc[:, [self.id_name,
                                               self.compound_name,
                                               self.retention_time_name]]

        except Exception as e:
            print("{} was not provided so it will continue without that information".format(e))
            self.name_data = self.raw_data.loc[:, [self.id_name,
                                                   self.compound_name]]

        self.name_data.drop_duplicates(subset=self.id_name, inplace=True)

        try:
            self.full_data = self.raw_data[[self.id_name,
                                            self.plant_name,
                                            self.target_variable_name,
                                            self.LOQ_flag_name]]

        except Exception as e:
            print("{} was not provided so it will continue without that information".format(e))
            self.full_data = self.raw_data.loc[:, [self.id_name,
                                                   self.plant_name,
                                                   self.target_variable_name]]

        return

    def curate_by_chemstruct(self, from_csv: bool = False, ignore_not_found: bool = False,
                             keep_stereochemistry: bool = False, verbose=False):
        """
        This method performs the following steps:

        1) Update names based on name dictionary (currently available for AMAR and AUS, else pass)
        2) Search for SMILES from Name using PubChemPy
        3) Manually add SMILES when necessary
        4) Canonicalize SMILES using rdkit
        5) Check for duplicates but this time based on canonical SMILES
        6) Create a combined ID and adds to full_data

        The method is applied to "name_data" which is the dataframe containing ID, Name and SMILES

        :param verbose: If verbose shows the curation steps.
        :param from_csv: If true, load existing csv files
        :param ignore_not_found: If True, compounds for which SMILES were not found using pubchempy
        are dropped from name_data for the rest of the workflow
        :param keep_stereochemistry: The default behaviour is to remove stereochemistry;
         if you want to keep it set to True
        :return: None
        """

        if from_csv:
            self.name_data = pd.read_csv(self.build_output_filename('smiles_info'), sep='\t')
            print("Existing file loaded from {}".format(self.name_data_tsv))
            return

        print("\n############# Curating based on chemical features  ############# ")

        # compound_name is replace with a list of names found in PubChem
        # The original names are kept as a column "RawNames"
        self.update_names(verbose=verbose)

        self.name_to_smiles(ignore_not_found=ignore_not_found, verbose=verbose)

        if self.tag == 'swe2_data':  # Temporary drop substances not found.
            self.name_data.dropna(subset=self.inchikey_name, inplace=True)

        # Make a copy of the original smiles to not overwrite them during curation
        self.name_data['Original_SMILES'] = self.name_data[self.smiles_name]

        if verbose:
            print("\n############# Check for disconnections in SMILES ############# ")
        # This  method takes a list and returns a list
        self.name_data[self.smiles_name] = Util.warn_disconnection(list(self.name_data[self.smiles_name]))

        if keep_stereochemistry:
            pass
        else:
            if verbose:
                print("\n############# Remove Stereochemistry ############# ")
            self.name_data[self.smiles_name] = self.get_rdkit_properties(remove_stereochemistry=True)

        if verbose:
            print("\n############# Canonicalizing SMILES ############# ")

        self.name_data['CanonicalSMILES'] = self.get_rdkit_properties(canonical_smiles=True)
        self.set_smiles_name('CanonicalSMILES')

        self.warn_duplicates(verbose=verbose)

        # Check point to save name_data with curated smiles
        self.name_data.to_csv(self.build_output_filename('smiles_info'), sep='\t', index=False)

        return

    def curate_by_target_variable(self, verbose=False, include_previous_criteria=False, only_above_LOQ=False):
        if self.full_data.empty:
            self.full_data = pd.read_csv(self.full_data_tsv, sep='\t')
            print('Full data was  loaded from {}'.format(self.full_data_tsv))

        print("\n############# Curating based on target variable values  ############# ")

        # This is to check the plants included in the analysis
        print("WWTPs included in analysis: {}".format(self.plant_list))
        self.full_data = self.full_data[self.full_data[self.plant_name].isin(self.plant_list)]

        if self.smiles_name not in self.full_data.columns:
            self.full_data = pd.merge(self.full_data, self.name_data[[self.smiles_name, self.id_name]], on=self.id_name)

        if only_above_LOQ:
            self.full_data = self.full_data[self.full_data[self.LOQ_flag_name] == 'False']
            # self.flag_LOQasEffluent(verbose=verbose)

        self.full_data['logB'] = np.log10(self.full_data.B)

        self.flag_not_enough_data(verbose=verbose)

        self.flag_suspected_formation(verbose=verbose)

        self.flag_stdev_level(verbose=verbose)

        if include_previous_criteria:
            self.filter_by_missing_data(verbose=verbose)  # changed to if NA in B ≥ 80% drop compound
            self.filter_by_logb()

        if self.tag != 'combined_data':
            self.full_data['dataset'] = self.tag
        else:
            self.full_data['dataset'] = 'multiple datasets'

        return

    def reduce_for_modelling(self, no_curation=True, only_enough_data=False,
                             no_formation=False, avoid_high_std=False, only_above_LOQ=False,
                             avoid_sorbing_and_volatile=False):
        """
        Reduces the dataset for modeling by applying various filtering criteria.

        This function sequentially applies the selected filtering criteria to the dataset to prepare it for modeling.
        It merges the full data with name data, filters based on the given parameters,
        and handles infinite and NaN values.

        :param no_curation: (bool) If True, no filtering will be applied, and the original dataset is used.
         Defaults to True.
        :param only_enough_data: (bool) If True, only data with sufficient samples (n >= 3) will be retained.
        :param no_formation: (bool) If True, substances suspected of formation (B > 1.2) will be removed.
        :param avoid_high_std: (bool) If True, substances with high standard deviation will be removed.
        :param only_above_LOQ: (bool) If True, only substances above the limit of quantification (LOQ) will be retained.
        :return: None. The function modifies the instance variable `self.cpd_data` to store the reduced dataset.
        """

        # Merge full data with name data on the specified ID
        self.full_data = pd.merge(self.full_data, self.name_data[[self.id_name,
                                                                  self.compound_name,
                                                                  'drop_duplicates']], on=self.id_name)


        # Print initial dataset information
        print("\n############# Selecting data for modeling  ############# ")
        print("Number of substances processed: {}".format(len(list(self.full_data[self.smiles_name].unique()))))
        self.full_data = self.full_data[self.full_data['drop_duplicates']]
        print("Number of substances after removing uncertain annotations: {}"
              .format(len(list(self.full_data[self.smiles_name].unique()))))

        # Initialize cpd_data with full_data
        self.cpd_data = self.full_data

        # Apply filtering criteria based on parameters
        if only_above_LOQ:
            self.curate_by_target_variable(only_above_LOQ=True)
            # choosing to use data above LOQ or not affects all other calculations, so it needs to be recalculated
            self.cpd_data = self.full_data
            print("Only above LOQ, current shape: {}".format(len(list(self.cpd_data[self.smiles_name].unique()))))

        if avoid_high_std:
            self.cpd_data = self.cpd_data[self.cpd_data['avoid_high_std']]
            print("Avoid substances with large standard deviation, current number of substances: {}"
                  .format(len(list(self.cpd_data[self.smiles_name].unique()))))

        if no_formation:
            self.cpd_data = self.cpd_data[self.cpd_data['not_suspected_formation']]
            print("Remove suspected formation, current number of substances: {}".
                  format(len(list(self.cpd_data[self.smiles_name].unique()))))

        if only_enough_data:
            self.cpd_data = self.cpd_data[self.cpd_data['enough_data']]
            print("Remove n<3, current shape: {}".format(len(list(self.cpd_data[self.smiles_name].unique()))))

        if avoid_sorbing_and_volatile:
            print("N substances before processing sorbing and volatile: {}".format(
                len(list(self.cpd_data[self.smiles_name].unique()))))
            self.flag_volatile_and_sorbing()
            # Keep only the ones confirmed to be 'not_volatile' & 'not_sorbind'
            self.cpd_data = self.cpd_data[self.cpd_data['not_volatile']==True]
            self.cpd_data = self.cpd_data[self.cpd_data['not_sorbing']==True]
            print("N substances after removing sorbing & volatile: {}".format(len(list(self.cpd_data[self.smiles_name].unique()))))


        elif no_curation:
            self.cpd_data = self.full_data
            print("No curation, data size: {}".format(len(list(self.cpd_data[self.smiles_name].unique()))))

        # Drop infinite values and rows with NaN values
        self.cpd_data.replace([np.inf, -np.inf], np.nan, inplace=True)  # Convert infinite values to NaN
        self.cpd_data.dropna(axis=0, inplace=True)  # Drop rows with NaN values

        # Print final dataset information
        print("Number of substances after dropping those with invalid values: {}"
              .format(len(list(self.cpd_data[self.smiles_name].unique()))))

        return

    def select_modeling_data(self, use_median=True, plant: str = False, dataset: str = False):
        if plant:
            self.cpd_data = self.cpd_data[self.cpd_data['plant'] == plant]
            print('plant {} was chosen for modeling'.format(plant))
            print(self.cpd_data.shape)
            return

        if dataset:
            self.cpd_data = self.full_data[self.full_data['dataset'] == dataset]
            print('dataset {} was chosen for modeling'.format(dataset))
            print(self.cpd_data.shape)

        if use_median:
            aggregation_functions = {
                self.target_variable_name: 'median',
                self.id_name: 'first',
                self.compound_name: 'first',
                'drop_duplicates': 'first',
                'enough_data': 'first',
                'not_suspected_formation': 'first',
                'avoid_high_std': 'first',
            }

            # Compute the standard deviation on the original data
            if self.target_variable_std_name != '':
                std_data = self.cpd_data.groupby(self.smiles_name)[self.target_variable_name].std().reset_index()
                std_data.rename(columns={self.target_variable_name: self.target_variable_std_name}, inplace=True)

            # Group by 'smiles_name', and apply the aggregation functions
            self.cpd_data = self.cpd_data.groupby(self.smiles_name).agg(aggregation_functions).reset_index()

            # Merge the computed standard deviation back to the aggregated DataFrame
            if self.target_variable_std_name != '':
                self.cpd_data = pd.merge(self.cpd_data, std_data, on=self.smiles_name, how='left')

            print('Median values were chosen for modeling, size : {}'.format(self.cpd_data.shape))

        else:
            print('all independent entries were chosen for modeling, size : {}'.format(self.cpd_data.shape))

        return

    def shape_aus_full_data(self):
        # Here I want to create name data and the IDs from the raw data
        self.name_data = self.raw_data[[self.compound_name, self.smiles_name, 'RT [min]']].copy()
        self.name_data[self.id_name] = self.name_data.index
        self.switch_to_alphanumeric_id(self.name_data)

        # Here I want to get full_data also from the raw data
        # The IDs should help me relate both name data and full data from the raw.

        dataframe = self.raw_data.copy()
        plant_df_dic = {}
#        default_plant_list = plants_group_2
        for plant_id in self.plant_list:
            plant_df = DataStructureWWTP.separate_data_by_plant(dataframe, str(plant_id))
            plant_df[self.id_name] = plant_df.index

            plant_df['plant'] = plant_id
            plant_df_dic[plant_id] = plant_df

        self.full_data = pd.concat(plant_df_dic, ignore_index=True)
        self.switch_to_alphanumeric_id(self.full_data)

        return

    @staticmethod
    def separate_data_by_plant(dataframe, plant_name):
        """
        Specific to AUS data set.
        Divides the "early" AUS full data in individual dataframes for each plant.
        It also drops compounds when influent values are missing.

        In AUS original raw file, influent effluent for each plant are provided in different columns.
        Each column is named as follows: "plant_name influent" OR "plant_name influent", for example "S2 influent"
        Thus the file is reorganized as independent dataframes with generic column names: 'influent' and 'effluent'.
        :param dataframe: A copy of AUS full data
        :param plant_name: The plant id that comes from datastructure.plant_list
        :return: dataframe with Names, SMILES, influent and effluent values
        """
        plant_df = copy.deepcopy(dataframe)
        influent = str(plant_name) + " influent"
        effluent = str(plant_name) + " effluent"
        plant_df.rename(columns={str(influent): 'influent', str(effluent): 'effluent'}, inplace=True)
        plant_df = plant_df[['influent', 'effluent']]
        return plant_df

    def format_amar_data(self):
        """
        Similarly to other 'format_dataset_data' methods the idea is to shape the raw file containing information
        of different plants into a single full_data dataframe. However, in the case of amar_data this curation
        was performed in a separated script because it is a very lengthy collection of rather repetitive functions
        See file new_amar_curation.ipynb
        """

        self.full_data = pd.read_csv(self.build_output_filename('raw_inf_eff'), sep='\t')
        self.name_data = self.raw_data[[self.compound_name, self.id_name, 'RT [min]']].copy()
        # Notice that aus_data would already have a self.smiles_name but not here.

        # self.name_data[self.id_name] = self.name_data.index
        # self.switch_to_alphanumeric_id(self.name_data)

    def format_snf_data(self):
        # self.name_data = self.raw_data
        inf_eff_df = pd.read_csv(self.build_output_filename('raw_inf_eff'), sep='\t')

        # Get the sampling week
        sampling_week_list = []
        for i in range(0, len(inf_eff_df.index)):
            sampling_week = inf_eff_df['Filename'][:][i][-1]
            sampling_week_list.append(sampling_week)
        inf_eff_df['sampling_week'] = sampling_week_list

        my_column_list = ['f.newcol_UchemID', 'Calculated.Amt', 'RT', 'LOQ',
                          'Calc.Amt_inf', 'LOQ_inf', 'WWTP', 'sampling_week', 'f.newcol_Names_input_Substance_list']
        rename_dict = {'f.newcol_UchemID': self.id_name, 'Calculated.Amt': 'effluent', 'LOQ': 'LOQ_eff',
                       'Calc.Amt_inf': 'influent', 'LOQ_inf': 'LOQ_inf', 'WWTP': self.plant_name,
                       'RT': self.retention_time_name, 'f.newcol_Names_input_Substance_list': self.compound_name}

        self.full_data = inf_eff_df[my_column_list].copy()
        self.full_data.rename(columns=rename_dict, inplace=True)

        # add alphanumeric id
        new_ID_list = []
        for index, row in self.full_data.iterrows():
            new_ID = str(self.tag) + '_' + str(row[self.id_name])
            new_ID_list.append(new_ID)
        self.full_data[self.id_name] = new_ID_list

        self.name_data = self.full_data[[self.compound_name, self.id_name, self.retention_time_name]].drop_duplicates()

        # I need to define name_data in terms of full_data
        # Afterwards, I want to drop all the information that does not belong to full_data
        self.full_data.drop(columns=[self.compound_name, self.retention_time_name], inplace=True)

        old_plant_label = list(self.full_data[self.plant_name])
        new_label_list = []
        for label in range(0, len(old_plant_label)):
            new_label = old_plant_label[label][1:-1]
            new_label_list.append(new_label)
        self.full_data[self.plant_name] = new_label_list

        return

    def format_galpro_data(self):
        self.full_data = pd.read_csv(self.build_output_filename('raw_inf_eff'), sep='\t')
        self.name_data = self.full_data[[self.compound_name, 'CompoundID']].drop_duplicates()
        self.name_data[self.id_name] = self.name_data.index

        # switch to alphanumerical ID
        new_ID_list = []
        for index, row in self.full_data.iterrows():
            new_ID = str(self.tag) + '_' + str(row['CompoundID'])
            new_ID_list.append(new_ID)
        self.full_data[self.id_name] = new_ID_list

        self.name_data = self.full_data[[self.compound_name, self.id_name,
                                         'CompoundID']].drop_duplicates()

    def filter_by_missing_data(self, verbose):
        """
        Normally not used (see 'flag_not_enough_data()' instead) but kept for documentation
        This function is also specific to aus_data and intends to exclude ,missing data
        As described by McLachlan(2022):
        'In order to generate a dataset suitable for exploring the influence of WWTP technology on breakthrough,
        two criteria for chemical inclusion were defined.
        First, breakthrough data were required for at least 80% (12) of the WWTPs
        to reduce the likelihood of data gaps biasing the results. '
        In the function '80% valid data points' is interpreted as 'number of missing data more than 20%'
        """

        remove_list = []
        self.name_data['missing_data'] = False
        for compound in list(self.name_data[self.compound_name]):
            if self.joint_data[self.joint_data[self.compound_name] == compound].loc[:, 'B(%)'].isna().sum() \
                    > int(self.joint_data[self.joint_data[self.compound_name] == compound].index.size)*0.2:
                remove_list.append(compound)

                # annotate name_data
                self.name_data['missing_data'] = np.where(self.name_data[self.compound_name] == compound, True,
                                                          self.name_data['missing_data'])
            else:
                pass
        if verbose:
            print("{} compounds will be removed due to insufficient data".format(len(remove_list)))

        return remove_list

    def filter_by_logb(self):
        """
        Normally not used (see 'flag_suspected_formation()' instead) but kept for documentation
        Checks the breakthrough values of each compound and keeps only the compounds for which log B is less than 2.2
        As described by McLachlan(2022):
        'Second, the 75th percentile of logB had to be <2.2 (corresponding to B < 158%)
        to exclude chemicals that were clearly formed during treatment.
        A total of 293 chemicals fulfilled these criteria'.
        """
        print("\n############# Annotating compounds suspected to be formed during treatment ############# ")

        suspected_formation = []
        self.name_data['suspected_formation'] = False
        i = 0
        for compound in list(self.name_data[self.compound_name]):
            if compound not in list(self.full_data[self.compound_name]):
                i += 1
                pass
            else:
                if np.nanpercentile(a=self.joint_data[self.joint_data[self.compound_name] == compound]['logB(%)'],
                                    q=75) >= 2.2:
                    print("{} is suspected to form during treatment".format(compound))
                    suspected_formation.append(compound)

                    # annotate name_data
                    self.name_data['suspected_formation'] = np.where(self.name_data[self.compound_name] == compound,
                                                                     True, self.name_data['suspected_formation'])
                else:
                    pass
        print("Target variable data is not yet available for {} compounds".format(i))
        return suspected_formation

    def flag_not_enough_data(self, na_threshold=3, verbose=False):
        """
        Flags compounds in the full_data dataframe that have insufficient data based on a given threshold.

        Args:
            na_threshold (int): The threshold for the number of data points required.
            verbose (bool): If True, prints the number of compounds that will be removed due to insufficient data.

        Returns:
            list: A list of compounds that will be removed due to insufficient data.
        """

        remove_list = []
        self.full_data['enough_data'] = True
        for compound in list(self.full_data[self.smiles_name].unique()):
            if self.full_data[self.full_data[self.smiles_name] == compound][self.target_variable_name].count() \
                    < na_threshold:
                remove_list.append(compound)

                # annotate name_data
                self.full_data['enough_data'] = np.where(self.full_data[self.smiles_name] == compound, False,
                                                         self.full_data['enough_data'])
            else:
                pass

        if verbose:
            print("{} compounds will be removed due to insufficient data".format(len(remove_list)))

        return remove_list

    def flag_suspected_formation(self, formation_threshold=1.20, verbose=False):
        """
        Flags compounds that are suspected to be formed during wastewater treatment, based on a given threshold.
        :param formation_threshold: Defined in terms of log10(B),
        The default is B = 1.20.
        :param verbose: If True returns the number of substances suspected to be formed.
        :return:
        """
        remove_list = []
        self.full_data['not_suspected_formation'] = True
        for compound in list(self.full_data[self.smiles_name].unique()):
            if np.nanmedian(self.full_data[self.full_data[self.smiles_name] == compound]
                            ['B']) > formation_threshold:
                remove_list.append(compound)

                # annotate full_data
                self.full_data['not_suspected_formation'] = np.where(self.full_data[self.smiles_name] == compound,
                                                                     False, self.full_data['not_suspected_formation'])
            else:
                pass

        if verbose:
            print("{} compounds will be removed due to suspected formation during treatment".format(len(remove_list)))

        return remove_list

    def flag_stdev_level(self, std_threshold=0.7, verbose=True):
        remove_list = []
        self.full_data['avoid_high_std'] = True
        for compound in list(self.full_data[self.smiles_name].unique()):
            if np.nanstd(self.full_data[self.full_data[self.smiles_name] == compound]
                         [self.target_variable_name]) > std_threshold:
                remove_list.append(compound)

                # annotate name_data
                self.full_data['avoid_high_std'] = np.where(self.full_data[self.smiles_name] == compound,
                                                            False, self.full_data['avoid_high_std'])
            else:
                pass

        if verbose:
            print("{} compounds have a high standard deviation".format(len(remove_list)))

        return remove_list

    def flag_volatile_and_sorbing(self, logHL_threshold=-5, logKoc_threshold=3.6, source='data'):

        # Load
        if source == 'data':
            self.opera_properties_tsv = '../data/wwtp-data/opera_properties_WWTP_combined_data.tsv'
        self.opera_properties = pd.read_csv(self.opera_properties_tsv, sep='\t')

        self.opera_properties.dropna(subset=['Opera_logHL', 'Opera_logKoc'], inplace=True)

        self.opera_properties['not_volatile'] = False
        self.opera_properties['not_sorbing'] = False

        # Flag (If condition is matched then True)
        # Should be interpreted as: if logHl < threshold then True, which means it is NOT volatile
        self.opera_properties['not_volatile'] = np.where(self.opera_properties['Opera_logHL'] < logHL_threshold,
                                                         True, False)

        # Should be interpreted as: if logKoc < threshold then True, which means it is NOT sorbing
        self.opera_properties['not_sorbing'] = np.where(self.opera_properties['Opera_logKoc'] < logKoc_threshold, True,
                                                        False)

        # Merge
        try:
            self.cpd_data = pd.merge(self.cpd_data, self.opera_properties, on=self.smiles_name, how='left')
        except Exception as e:
            print('Merging the opera_properties frame did not work: {}'.format(e))

        return

    def flag_LOQasEffluent(self):
        """ Flag those values for which the LOQ was used to calculate effluent concentrations.
            Cases where the influent is below the LOQ are never used.
        """
        # Missing values of B are flagged before being replaced by upper limit.
        if self.tag == 'swe2_data':
            self.full_data['LOQasEffluent'] = self.full_data.B.isna()

        elif self.tag == 'amar_data':
            pass

        else:
            print('It has not being defined how to ')
            pass

        return


    def get_my_names_dict(self):
        filename = '../config/wwtp_names.yml'
        with open(filename, 'r') as file:
            wwtp_names_dict = yaml.safe_load(file)
        if self.tag == 'aus_data':
            return wwtp_names_dict.aus_names_dict

        if self.tag == 'amar_data':
            from pepper.wwtp_names import amar_names_dict
            return amar_names_dict
        if self.tag == 'snf_data':
            from pepper.wwtp_names import snf_names_dict
            return snf_names_dict

        if self.tag == 'aus_amar_data':
            from pepper.wwtp_names import aus_amar_names_dict
            return aus_amar_names_dict

        if self.tag == 'galpro_data':
            from pepper.wwtp_names import galpro_names_dict
            return galpro_names_dict

        if self.tag == 'swe2_data':
            print('For now I am updating with AMAR_dict')
            from pepper.wwtp_names import amar_names_dict
            return amar_names_dict

        else:
            print('names_dic not available')
            return

    def update_names(self, verbose: bool = False):
        """
        Updates names as provided in raw data files to names as found in PubChem.
        The names are updated based on a manually compiled dictionary.
        :param verbose: If true, prints the names that were updated.
        :return:
        """

        name_list = self.name_data[self.compound_name]

        #
        my_names_dic = self.get_my_names_dict()
        if not my_names_dic:
            print('Names were not updated because the dictionary of new names has not been provided ')
            return

        pubchem_names_list = []
        i = 0
        for name in name_list:
            if name in my_names_dic.keys():
                i += 1
                if verbose:
                    print(name)
                    print("found it! Updated to:")
                    print(my_names_dic.get(name))
                pubchem_names_list.append(my_names_dic.get(name))
            else:
                pubchem_names_list.append(name)

        if verbose:
            print("\n############# Updating names to match PubChem ############# ")
            print("{} names were updated".format(i))

        # Update name to keep raw names and not overwrite them
        self.name_data.rename(columns={self.compound_name: 'RawNames'}, inplace=True)

        # Save the new names as the new self.compound_name column
        self.name_data[self.compound_name] = pubchem_names_list

        return

    def manual_fix_smiles(self):
        """
        Kept only for documentation. This compound was not found in PubChem at the time
        Manual annotation of compounds of the AUS dataset  which were not found on the PubChem database """

        # compound missing in pubchem
        cpd_not_found = '(2E)-2-[(2,1,3-benzothiadiazol-4-ylamino)methylidene]-4,4-dimethyl-3-oxopentanenitrile'
        correct_smiles = 'CC(C)(C)C(=O)C(=CNC1=CC=CC2=NSN=C21)C#N'
        self.name_data[self.smiles_name] = np.where(self.name_data[self.compound_name] == cpd_not_found,
                                                    correct_smiles, self.name_data[self.smiles_name])

        return

    @staticmethod
    def get_k(breakthrough, k_hyd):
        return k_hyd * ((1-breakthrough)/breakthrough)

    def add_k_bio(self, k_hyd):  # todo: Check if this is really going to be implemented as some point; else remove
        """

        """
        print("Calculating kinetic constants from breakthrough assuming steady state")

        if 'B' in self.full_data.columns:
            self.full_data['k'] = DataStructureWWTP.get_k(self.full_data.B, k_hyd)
        elif 'B(%)' in self.full_data.columns:
            B = self.full_data['B(%)'] / 100
            self.full_data['k'] = DataStructureWWTP.get_k(B, k_hyd)
        else:
            print("k was not calculated. B or B (%) must be calculated before getting k")

    def warn_duplicates(self, verbose=True):
        """
        Get unique counts of compounds by index, Name & SMILES individually.
        Warn if there are duplicates and if there are, create a dataframe of duplicates.
        Finally, flag duplicates in the name_data dataframe (see method flag_duplicates()).
        :param verbose: If true, prints the names that were updated.
        :return:
        """

        # Check counts in each plant

        unique_names = self.name_data[self.compound_name].unique().size
        unique_SMILES = self.name_data[self.smiles_name].unique().size

        if verbose:
            print("\n############# Check duplicates ############# ")
            print("Number of unique names: {}".format(unique_names))
            print("Number of unique substances based on canonicalized SMILES strings: {}".format(unique_SMILES))
            if unique_names != unique_SMILES:
                print("\n############# WARNING there seems to be duplicates ############# ")

        self.duplicates_df = self.name_data[self.name_data.duplicated(subset=self.smiles_name, keep=False).copy()]

        if self.duplicates_df.empty:
            self.flag_duplicates(ignore=True)
        else:
            print(list(self.duplicates_df[self.compound_name]))
            # Flag duplicates requires retention time data which is not always available
            if self.retention_time_name in self.duplicates_df.columns:
                self.flag_duplicates(verbose=verbose, ignore=False)
            else:
                print('There is no information about retention time so duplicates will not be processed.'
                      'We suggest checking "self.duplicates_df" to make a decision and reprocess the data.')
                pass

        return

    def flag_duplicates(self, ignore=False, verbose=False):
        """
        Check the retention time to decide whether the peak assignment could have been wrong.
        IF the retention time differs in more than one percent, drop both.
        :param ignore: if true, do not flag duplicates but create columns for drop_duplicates and Combined_ID
        for consistency.
        :param verbose:
        :return:
        """
        if ignore:
            self.name_data['drop_duplicates'] = True
            try:
                self.name_data['Combined_ID'] = self.name_data[self.id_name]
            except Exception as e:
                print("({}. Combined_IDs being ignored; Check reasons".format(e))
            return
        else:
            duplicates_RT_df = self.split_duplicates()
            duplicates_RT_df['RT_diff'] = (np.abs(duplicates_RT_df['RT [min]_1']-duplicates_RT_df['RT [min]_2']) /
                                           duplicates_RT_df['RT [min]_1'])*100

            duplicates_drop_list = list(duplicates_RT_df[duplicates_RT_df['RT_diff'] > 1][self.id_name+'_1'])+list(
                duplicates_RT_df[duplicates_RT_df['RT_diff'] > 1][self.id_name+'_2'])

            drop_list = []
            for index, row in self.name_data.iterrows():
                if row.ID in duplicates_drop_list:
                    a = False
                else:
                    a = True
                drop_list.append(a)
            self.name_data['drop_duplicates'] = drop_list

            ref_ID_df = self.create_ref_ID(duplicates_RT_df)
            if verbose:
                print(ref_ID_df)
            self.name_data = pd.merge(self.name_data, ref_ID_df, how='outer', on=self.smiles_name)

            # when there is no combined_ID because there are no duplicates, use the ID
            self.name_data.Combined_ID.fillna(self.name_data.ID, inplace=True)

            return

    def merge_duplicates(self, verbose=False):
        # I prefer to merge on join data because name_data contains the relevant information
        # for handling the duplicates correctly, however, this information is not
        # reflected in full data, yet.
        # Thus, I merge then I handle the duplicates and then I redefine name_data and full_data

        self.joint_data = pd.merge(self.name_data, self.full_data, on='ID')

        aggregation_functions = {
            'B': 'mean',
            # 'ID': 'first',
            # 'Original_SMILES': 'first',
            # 'RT [min]': 'first',
            # 'InChIKey': 'first',
            'dataset': 'first',
            'LOQasEffluent': 'first',
            # 'ClogP': 'first',
            'drop_duplicates': 'first',
            self.compound_name: 'first',
            self.smiles_name: 'first'}

        self.joint_data = self.joint_data.groupby(['Combined_ID', self.plant_name]).agg(aggregation_functions).reset_index()

        # Group by names, IDs and the drop_duplicates tag
        self.name_data = self.joint_data[['Combined_ID', self.compound_name, self.smiles_name, 'drop_duplicates']]
        # This definition gets the right info but expanded to all plants.
        # Now that I have dealt with duplicates, if I keep only the unique entries it
        # will have the right size (one entry per molecule).
        self.name_data.drop_duplicates(inplace=True)
        if verbose:
            print(self.name_data.shape)

        self.full_data = self.joint_data[['Combined_ID', 'plant', 'B', 'dataset', 'LOQasEffluent']]
        return

    def create_ref_ID(self, dataframe):
        ref_id_df = dataframe[[self.id_name+'_1', self.id_name+'_2', self.smiles_name]].copy()
        id_ref = []
        for index, row in ref_id_df.iterrows():
            a = str(str(row.ID_1) + '+' + str(row.ID_2))
            id_ref.append(a)
        ref_id_df['Combined_ID'] = id_ref
        return ref_id_df[['Combined_ID', self.smiles_name]]

    def split_duplicates(self):
        smiles_list = []
        row_1 = []
        row_2 = []
        for index, row in self.duplicates_df[[self.id_name, self.smiles_name, self.retention_time_name]].iterrows():
            RT = row
            a = row[self.smiles_name]
            if a not in smiles_list:
                row_1.append(RT)
                smiles_list.append(a)
            else:
                row_2.append(RT)
        return pd.merge(pd.DataFrame(row_2), pd.DataFrame(row_1), on=self.smiles_name, suffixes=('_1', '_2'))

# Static methods to study the kinetics of each compound considering breakthrough values
# single code to get all figures

    @staticmethod
    def get_k_steady_state(breakthrough, time_value):
        return time_value * ((1 - breakthrough) / breakthrough)

    @staticmethod
    def get_k_single_compartment(breakthrough, time_value):
        return -np.log(breakthrough) / time_value

    # get B
    # Steady state
    @staticmethod
    def get_B_from_k(k, time_value):
        return time_value / (k + time_value)

    # Single compartment
    @staticmethod
    def get_B_from_k_single(k, time_value):
        return np.exp(-k * time_value)

    @staticmethod
    def get_k_profile(x_range, time_value=0.1, kinetics_type='steady_state',
                      get_B=True, logx=False, logy=False):
        k = np.linspace(x_range[0], x_range[1], 10000000)
        B = np.linspace(x_range[0], x_range[1], 10000000)

        if get_B:
            y_label = 'breakthrough'
            if kinetics_type == 'steady_state':
                plt.plot(k, DataStructureWWTP.get_B_from_k(k, time_value), color='red')
            elif kinetics_type == 'single_compartment':
                plt.plot(k, DataStructureWWTP.get_B_from_k_single(k, time_value), color='red')
            else:
                print("kinetics_type {} is not supported".format(kinetics_type))

        else:
            y_label = 'kinetic constant'
            if kinetics_type == 'steady_state':
                plt.plot(B, DataStructureWWTP.get_k_steady_state(B, time_value), color='red')
            elif kinetics_type == 'single_compartment':
                plt.plot(k, DataStructureWWTP.get_k_single_compartment(B, time_value), color='red')

        if logx:
            plt.xscale('log')
        if logy:
            plt.yscale('log')

        plt.ylabel(str(y_label))
        plt.show()
        return plt

    def get_canonicalized_smiles(self):
        # todo: overrides method in DataStructure.
        #  The only reason is that name_data currently does not exist in DataStructure

        """
        Util.canonicalize_smiles canonicalizes SMILES strings one by one.
        This method creates a list of canonical SMILES by calling the Util method to each element on full_data.
        :return: List of canonical SMILES
        """
        new = []
        for index, row in self.name_data.iterrows():
            full_smiles = row[self.smiles_name]
            canonical_smiles = Util.canonicalize_smiles(full_smiles)
            new.append(canonical_smiles)
        return new

    def get_rdkit_properties(self, canonical_smiles=False,
                             clogp=False,
                             remove_stereochemistry=False):
        """
        Util.canonicalize_smiles canonicalizes SMILES strings one by one.
        This method creates a list of canonical SMILES by calling the Util method to each element on name_data.
        :return: List of canonical SMILES
        """
        new = []
        for index, row in self.name_data.iterrows():
            full_smiles = row[self.smiles_name]
            if canonical_smiles:
                canonical_smiles = Util.canonicalize_smiles(full_smiles)
                new.append(canonical_smiles)
            elif clogp:
                calculated_logp = Util.get_clogp(full_smiles)
                new.append(calculated_logp)
            elif remove_stereochemistry:
                simple_smiles = Util.remove_stereochemistry(full_smiles)
                new.append(simple_smiles)
            else:
                print("Select property")

        return new

    def correct_by_LODQ(self, verbose=False):
        """
        This method checks the LOQ and LOD values for both influent and effluent.
        If the influent is less than the LOQ the influent is set to "np.nan" and so no B values will be calculated
        If the effluent is below the LOQ but above the LOD then the LOD value is kept and a flag is added.
        If the effluent is below the LOD then the LOD is used and another flag is added.
        These rules were selected for the AMAR dataset.
        In the case of SNF dataset compounds are flagged if effluent is below the LOQ.
        In the case of the AUS dataset, only values above the LOQ are reported.
        In the case of Yijing data the flag is provided.
        """

        self.full_data[self.LOQ_flag_name] = False
        if verbose:
            print("\n############# Correcting breakthroughs below LOQ ############# ")
        self.full_data[self.LOQ_flag_name] = np.where(self.full_data['effluent'] < self.full_data['LOQ'],
                                                      '<LOQ', self.full_data[self.LOQ_flag_name])
        self.full_data[self.LOQ_flag_name] = np.where(self.full_data['effluent'] < self.full_data['LOD'],
                                                      '<LOD', self.full_data[self.LOQ_flag_name])

        self.full_data['effluent'] = np.where(self.full_data['effluent'] < self.full_data['LOD'],
                                              self.full_data['LOD'], self.full_data['effluent'])

        self.full_data[self.LOQ_flag_name] = np.where(self.full_data['influent'] < self.full_data['LOQ'],
                                                      'invalid influent', self.full_data[self.LOQ_flag_name])

        self.full_data['influent'] = np.where(self.full_data['influent'] < self.full_data['LOQ'],
                                              np.nan, self.full_data['influent'])

        return

    def snf_correct_by_LODQ(self, verbose=False):
        self.full_data[self.LOQ_flag_name] = False
        if verbose:
            print("\n############# Correcting breakthroughs below LOQ ############# ")

        # more or less: if effluent is below the LOQ use the LOQ and flag it, else keep the effluent value.
        self.full_data[self.LOQ_flag_name] = np.where(self.full_data['effluent'] < self.full_data['LOQ_eff'],
                                                      '<LOQ', self.full_data[self.LOQ_flag_name])

        self.full_data['effluent'] = np.where(self.full_data['effluent'] < self.full_data['LOQ_eff'],
                                              self.full_data['LOQ_eff'], self.full_data['effluent'])

        # If the influent is below the LOQ set as np.nan to drop later, else keep the influent value.
        # We don't want to consider data with influent below the LOQ because we consider it too uncertain.
        self.full_data['influent'] = np.where(self.full_data['influent'] < self.full_data['LOQ_inf'],
                                              np.nan, self.full_data['influent'])
        return

    def merge_snf_entries(self):
        # todo: I would prepare to calculate the breakthrough first and then get the mean of that breakthrough
        #  instead of getting the breakthrough from mean influent and effluent
        self.full_data['B'] = self.full_data.effluent/self.full_data.influent

        aggregation_functions = {
            'B': 'median',
            'LOQ_inf': 'first',
            'LOQ_eff': 'first',
            self.LOQ_flag_name: 'first'}

        # Get a single, mean entry per triplicated
        self.full_data = self.full_data.groupby(['plant', self.id_name]).agg(aggregation_functions).reset_index()

    def switch_to_alphanumeric_id(self, dataframe):
        new_ID_list = []
        for index, row in dataframe.iterrows():
            new_ID = self.tag + '_' + str(row[self.id_name])
            new_ID_list.append(new_ID)
        dataframe[self.id_name] = new_ID_list
