from copy import deepcopy

import pandas as pd
from pepper_lab.pepper import Pepper
from pepper_lab.metadata import *
from pepper_lab.util import *
from pepper_lab.visualize import Visualize

from enviPath_python import enviPath
from enviPath_python.objects import *

import getpass
from tqdm import tqdm
from sklearn.metrics import mean_squared_error, r2_score


class DataStructure(Pepper):
    def __init__(self, pep: Pepper):
        super().__init__()

        # General attributes for all datastructure classes
        self.set_data_directory(os.path.join(pep.data_directory, 'data_structure'))
        # original data as loaded from source
        self.raw_data = pd.DataFrame()
        self.raw_data_tsv = self.build_output_filename('raw_data')
        # Clean data that identifies each unique substance in the dataset
        self.name_data = pd.DataFrame()
        self.name_data_tsv = self.build_output_filename('name_data')

        # current, annotated data
        self.full_data = pd.DataFrame()
        self.full_data_tsv = self.build_output_filename('full_data')
        # reduced data with extended information
        self.cpd_data = pd.DataFrame()
        self.cpd_data_tsv = self.build_output_filename('cpd_data')
        # data used for modelling with only with standardized column names
        self.model_data = pd.DataFrame()
        self.model_data_tsv = self.build_output_filename('model_data_' + self.setup_name + '_' + self.curation_type)
        # Dict to store data
        self.data_dict = {}

        # Attributes that we want to keep from the pepper object
        self.pepper = pep
        self.tag = pep.get_tag()
        self.data_type = pep.get_data_type()
        self.compound_name = pep.get_compound_name()
        self.smiles_name = pep.get_smiles_name()
        self.target_variable_name = pep.get_target_variable_name()
        self.target_variable_std_name = pep.get_target_variable_std_name()
        self.id_name = pep.get_id_name()
        self.random_state = pep.get_random_state()
        self.plant_name = 'plant'
        self.target_variable_list = ['endpoint']
        self.inchikey_name = 'InChIKey'

        # enviPATH
        self.instance_host = "https://envipath.org"
        self.envipath_package = ''

    def set_envipath_package(self, uri):
        self.envipath_package = uri

    def set_target_variable_list(self, target_variable_list):
        """
        Specify the list of target variables that will be considered. 
        This will allow to calculate values for all these variables.
        However, the model will be trained and tested according to a single variable specified in target_variable_name.
        :param target_variable_list:
        """
        self.target_variable_list = target_variable_list

    def set_curation_type(self, curation_type):
        self.curation_type = curation_type
        # update the model_data_tsv when updating the curation type
        self.model_data_tsv = self.build_output_filename('model_data_' + self.setup_name + '_' + self.curation_type)

    def get_target_variable(self):
        return self.model_data[self.target_variable_name]

    def get_smiles(self):
        return self.full_data[self.smiles_name]

    def load_raw_data(self, source='pepper_data'):
        """
        Loads the raw data from pepper_data, data or enviPath. todo: remove and use load_data instead
        :param source: possible values: "pepper_data" (existing data file),
        "data" (from raw_data folder), "enviPath" (download from envipath.org)
        """
        print("\n############# Loading raw data ############# ")

        if source == "pepper_data": # todo : change naming
            assert os.path.exists(self.raw_data_tsv), "Error: file {} does not exist".format(self.raw_data_tsv)
            self.raw_data = pd.read_csv(self.raw_data_tsv, sep='\t')
        elif source == "data":
            file_string = ('raw_data' + '_{}_{}' + '.tsv').format(self.data_type, self.tag)
            self.raw_data = pd.read_csv(os.path.join('..', 'data', file_string), sep='\t')
        elif source == "enviPath":
            self.load_raw_data_from_enviPath()
        else:
            raise Exception("Unknown source")

    def load_data(self, data_type, source=False):
        """
        Loads the raw data from pepper_data, data or enviPath.
        :param data_type: possible values: 'raw_data', 'full_data', 'name_data', 'cpd_data', 'model_data'
        :param source: Default (None) is to look for the local directory; use data when using a virtual environment
        """
        print("\n############# Loading data ############# ")

        data_tsv = getattr(self, '{}_tsv'.format(data_type))

        if source == "data":
            tsv_name = (data_type + '_' + self.setup_name + '_' + self.data_type + '_' + self.tag + '.tsv')
            data_tsv = os.path.join('..', 'data', 'wwtp-data', tsv_name)


        assert os.path.exists(data_tsv), "Error: file {} does not exist".format(data_tsv)
        my_data = pd.read_csv(data_tsv, sep='\t')

        if my_data is not None:
            setattr(self, data_type, my_data)
            print('{} loaded from {}'.format(data_type, data_tsv))
        else:
            print("Data could not be loaded from {}".format(data_tsv))


    def get_model_data(self):
        return self.model_data

    def get_name_data(self):
        return self.name_data

    @staticmethod
    def check_for_kinetics(addinfo):
        try:
            addinfo.get_halflife().get_value()
        except AttributeError:
            return False
        else:
            return True

    # functions
    def load_raw_data_from_enviPath(self):

        eP = enviPath(self.instance_host)
        pkg = Package(eP.requester, id=self.envipath_package)
        pathways = pkg.get_pathways()
        for path in tqdm(pathways[:3]):
            print(pathways.index(path), path.get_id())
            for node in path.get_nodes():
                scenarios = node.get_scenarios()
                for scenario in scenarios:
                    # print(scenario.get_id())
                    full_scenario = Scenario(eP.requester, id=scenario.get_id())
                    temp_add_info = full_scenario.get_additional_information()
                    add_info = {ai.name: ai for ai in temp_add_info}
                    description = full_scenario.get_description()  # to obtain high and low OC information
                    if any([True if (isinstance(obj, RateConstantAdditionalInformation) or
                                     isinstance(obj, HalfLifeAdditionalInformation)) else False for obj in add_info.values()]):
                        # load all necessary data form enviPath
                        compound = CompoundStructure(eP.requester, id=node.get_default_structure().get_id())
                        metadata = Metadata(add_info, description)
                        try:
                            spike_smiles = CompoundStructure(
                                eP.requester, id=add_info.get_spikecompound().get_compoundLink()).get_smiles()
                        except:
                            spike_smiles = ''

                        self.data_dict = metadata.get_scenario_information(self.data_dict, scenario,
                                                                           compound, self.data_type,
                                                                           spike_smiles, description)

        # save data
        self.raw_data = pd.DataFrame(self.data_dict)
        self.raw_data.to_csv(self.raw_data_tsv, sep='\t', index=False)

    def curate_smiles(self, from_csv=False):
        print("\n############# SMILES curation ############# ")
        if from_csv:
            self.full_data = pd.read_csv(self.full_data_tsv, sep='\t')
            print("Existing file loaded from {}".format(self.full_data_tsv))
            return
        cropped_smiles, is_composite = self.get_cropped_smiles()
        self.full_data['cropped_SMILES'] = cropped_smiles
        self.full_data['is_composite'] = is_composite
        self.full_data['canonical_SMILES'] = self.get_canonicalized_smiles()
        self.full_data['cropped_canonical_SMILES'] = self.get_cropped_canonicalize_smiles()
        self.full_data['cropped_canonical_SMILES_no_stereo'] = self.get_cropped_canonicalize_smiles_no_stereo()

        self.full_data.rename(columns={self.smiles_name: 'original_SMILES'})  # set smiles to work with
        self.full_data[self.smiles_name] = self.full_data['cropped_canonical_SMILES_no_stereo']  

        self.full_data.to_csv(self.full_data_tsv, sep='\t', index=False)

    def create_modelling_input(self, from_csv=False, no_curation=False, include_plant=False):
        """
        Creates a dataframe with the necessary data to create the models.
        That is, IDs, SMILES strings and names of target compounds and values of the target variable.
        It also shuffles the data.
        :param from_csv:

        Parameters
        ----------
        :param no_curation: bool default False If true, no curation is performed
        :param include_plant: bool default False, use False when using median for modeling.
        If instead one would like to use the individual values for each compound in each plant,
        then set include_plant to True.
        """
        if no_curation:
            self.model_data = self.cpd_data
            return

        if from_csv:
            assert os.path.exists(self.model_data_tsv), "Error: file {} does not exist".format(self.model_data_tsv)
            self.model_data = pd.read_csv(self.model_data_tsv, sep='\t')
            return

        self.cpd_data = self.cpd_data.sample(frac=1, random_state=42, ignore_index=True)

        self.model_data = self.cpd_data.loc[:, [self.id_name,
                                                self.smiles_name,
                                                self.compound_name,
                                                self.target_variable_name]]

        # Add std to target variable if available
        if self.target_variable_std_name != '':
            self.model_data[self.target_variable_std_name] = self.cpd_data[self.target_variable_std_name]

        # Add plant identifier if necessary
        if include_plant:
            self.model_data[self.plant_name] = self.cpd_data[self.plant_name]

        self.model_data.to_csv(self.model_data_tsv, sep='\t', index=False)

    def randomize_y(self):
        random_y = (self.model_data[self.target_variable_name].sample(frac=1, random_state=42,
                                                                      replace=False, ignore_index=False)).values

        self.model_data[self.target_variable_name] = random_y
        
    def get_cropped_smiles(self):
        new_composite = []
        new_smiles = []
        full_smiles = ''
        cropped_smiles = ''
        is_composite = False
        for index, row in self.full_data.iterrows():
            if row[self.smiles_name] != full_smiles:
                full_smiles = row[self.smiles_name]
                cropped_smiles = full_smiles
                is_composite = False
                if '.' in full_smiles:
                    smiles_list = full_smiles.split('.')
                    smiles_list.sort(key=len, reverse=True)
                    cropped_smiles = smiles_list[0]
                    # special case of fluroxypyr ester + acid: keep acid part
                    if full_smiles == 'CCCCCCC(C)OC(=O)COc1nc(F)c(Cl)c(N)c1Cl.Nc1c(Cl)c(F)nc(OCC(=O)O)c1Cl':
                        cropped_smiles = 'Nc1c(Cl)c(F)nc(OCC(=O)O)c1Cl'
                    is_composite = True
            new_smiles.append(cropped_smiles)
            new_composite.append(is_composite)
        return new_smiles, new_composite

    def get_canonicalized_smiles(self):
        """
        Util.canonicalize_smiles canonicalizes SMILES strings one by one.
        This method creates a list of canonical SMILES by calling the Util method to each element on full_data.
        :return: List of canonical SMILES
        """
        new = []
        for index, row in self.full_data.iterrows():
            full_smiles = row[self.smiles_name]
            canonical_smiles = Util.canonicalize_smiles(full_smiles)
            new.append(canonical_smiles)
        return new

    def get_cropped_canonicalize_smiles(self):
        new = []
        cropped_smiles = ''
        cropped_canonical_smiles = ''
        for index, row in self.full_data.iterrows():
            if row['cropped_SMILES'] != cropped_smiles:
                cropped_smiles = row['cropped_SMILES']
                cropped_canonical_smiles = Util.canonicalize_smiles(cropped_smiles)
            new.append(cropped_canonical_smiles)
        return new

    def get_cropped_canonicalize_smiles_no_stereo(self):
        new = []
        cropped_canonical_smiles = ''
        cropped_canonical_smiles_no_stereo = ''
        for index, row in self.full_data.iterrows():
            if row['cropped_canonical_SMILES'] != cropped_canonical_smiles:
                cropped_canonical_smiles = row['cropped_canonical_SMILES']
                cropped_canonical_smiles_no_stereo = Util.canonicalize_smiles(
                    Util.remove_stereo_info(cropped_canonical_smiles))
            new.append(cropped_canonical_smiles_no_stereo)
        return new

    def fetch_few_compounds(self, n_compounds):
        self.model_data = self.model_data.iloc[:n_compounds]


    def experimental_performance_simulation(self,
                                            number_of_samples: int = 10, reported_experimental_value: str = 'DT50_log'):
        """
        @param type:  'experimental_values' or 'samples_from_bayesian_distribution'
        @param number_of_samples: number of samples for performance calculation
        @param reported_experimental_value: for typ='experimental_values', the column where the reported target variables are reported
        """
        print("--> Run performance simulation. Number of samples: {}".format(number_of_samples))
        r2_exp, rmse_exp = self.get_performance_from_experimental_values(number_of_samples, reported_experimental_value)
        r2_dist, rmse_dist = self.get_performance_from_distribution(number_of_samples)
        df = pd.DataFrame()
        df['R2_exp'] = r2_exp
        df['RMSE_exp'] = rmse_exp
        df['R2_dist'] = r2_dist
        df['RMSE_dist'] = rmse_dist
        print(df.describe())
        v = Visualize(self, "experimental_performance_simulation")
        v.plot_experimental_performance_simulation(df)


    def get_performance_from_distribution(self, number_of_samples):
        samples = []
        true_mean = self.cpd_data[self.target_variable_name]
        np.random.seed(self.pepper.random_state)
        for index, row in self.cpd_data.iterrows():
            mean = row[self.target_variable_name]
            std = row[self.target_variable_std_name]
            samples.append(np.random.normal(loc=mean, scale=std, size=number_of_samples))
        r2 = []
        rmse = []
        for index, sample in enumerate(np.array(samples).T): # iterate through 100 samples
            # split sample into 5 and determine performance on each subset.
            sample_index = np.arange(len(samples))
            indices_list = Util.split_function(sample_index, 5, seed=index)
            for indices in indices_list:
                rmse.append(mean_squared_error(true_mean[indices], sample[indices]))
                r2.append(r2_score(true_mean[indices], sample[indices]))
        return r2, rmse

    def get_performance_from_experimental_values(self, number_of_samples, reported_experimental_value):
        samples = []
        true_mean = []
        np.random.seed(self.pepper.random_state)
        for index, row in self.cpd_data.iterrows():
            this_data = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            experimental_values = this_data[reported_experimental_value].values
            if len(experimental_values) >= 3:
                samples.append(np.random.choice(list(experimental_values), size=number_of_samples))
                true_mean.append(row[self.target_variable_name].real)
        print("Found {} compounds with 3 or more data points to be considered for analysis (out of {}).".format(len(samples), len(self.cpd_data)))
        r2 = []
        rmse = []
        true_mean = np.array(true_mean)
        for index, sample in enumerate(np.array(samples).T):
            # determine size of subsets for performance calculation. Should be in line with test set size in 5-fold CV
            subset_length = round(len(self.cpd_data) / 5)
            for subset_id in np.arange(5):
                sample_index = np.arange(len(samples))
                indices = np.random.choice(sample_index, size=subset_length)
                rmse.append(mean_squared_error(true_mean[indices], sample[indices]))
                r2.append(r2_score(true_mean[indices], sample[indices]))
        return r2, rmse

    def analyze_parameter_distributions(self, reported_endpoint_name):
        """
        Print statistics and plot distributions of experimental and environmental parameters
        """
        print("\n############# Analyze experimental parameter distribution ############# ")
        header_list = deepcopy(self.experimental_parameter_names)
        header_list.append(reported_endpoint_name)
        df_params = self.full_data.loc[:, header_list]
        # Save statistics
        df_params.describe()
        # Visualize distributions
        v = Visualize(self, 'analyze_distributions')
        v.plot_experimental_parameter_distribution(df_params, reported_endpoint_name)

