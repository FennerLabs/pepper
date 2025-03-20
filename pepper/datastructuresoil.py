import os
from pepper.pepper import Pepper
# from pepper.metadata import *
from pepper.util import *
from pepper.datastructure import DataStructure
from pepper.bayesian import *
from pepper.visualize import Visualize


class DataStructureSoil(DataStructure):
    def __init__(self, pep: Pepper):
        super().__init__(pep)
        self.set_data_directory(os.path.join(pep.data_directory,'data_structure','soil'))
        self.smiles_name = pep.get_smiles_name()
        self.target_variable_name = pep.get_target_variable_name()
        self.target_variable_std_name = pep.get_target_variable_std_name()
        self.compound_name = pep.get_compound_name()
        self.id_name = pep.get_id_name()

        self.envipath_package = 'https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456'
        self.data_type = 'soil'
        self.data_dict = {}

        # These are general, yet, they must be initialized with the info of the new object
        # Can we avoid defining all these for every subclass?
        self.raw_data_tsv = self.build_output_filename('raw_data')
        self.full_data_tsv = self.build_output_filename('full_data')
        self.cpd_data_tsv = self.build_output_filename('cpd_data')
        self.model_data_tsv = self.build_output_filename('model_data')
        self.cpd_data_description_file = self.build_output_filename('cpd_data_description')

        # Soil specific attributes
        self.spike_compound_dictionary = {}


    def curate_annotate(self, from_csv: bool = False, from_paper: bool = False):
        """
        Curate SMILES, log-transform values, calculate bayesian-inferred mean/std of target variable,
        and curate information on half-lives.
        :param from_csv: If true, load existing csv files
        :return:
        """
        if from_csv:
            if from_paper:
                self.full_data = pd.read_csv(os.path.join('..', 'data', 'soil', 'full_data_soil_all_data.tsv'), sep='\t')
                return
            else:
                self.full_data = pd.read_csv(self.full_data_tsv, sep='\t')
                print("Existing curated file loaded from {}".format(self.full_data_tsv))
                return

        self.full_data = self.raw_data
        self.curate_smiles()

        self.curate_halflife_data()

    def transform_values(self, from_csv = False):
        print("\n############# Value transformation - log transformation and bayesian inference ############# ")
        if from_csv:
            self.full_data = pd.read_csv(self.full_data_tsv, sep='\t')
            print("Existing file loaded from {}".format(self.full_data_tsv))
            return

        # Index compounds by smiles identity
        self.full_data[self.id_name] = self.index_compounds()
        self.full_data.sort_values(by=self.id_name)

        self.full_data['DT50_log'] = Util.log_transform(self.full_data['reported_DT50'])
        self.full_data['soil_texture_main'] = self.get_main_soil_texture()
        self.full_data['CEC_log'] = Util.log_transform(self.full_data['CEC'])
        self.full_data['biomass_log'] = Util.log_transform(self.full_data['biomass'])
        self.full_data['OC_log'] = Util.log_transform(self.full_data['OC'])
        self.full_data['DT50_count'] = self.count_halflives()
        self.full_data['DT50_log_median'] = self.get_hl_median()  # median
        self.full_data['DT50_gmean'] = self.get_geometric_mean()  # geometric mean
        self.full_data['DT50_log_gmean'] = Util.log_transform(self.get_geometric_mean())  # log of geometric mean
        self.full_data['DT50_log_std'] = self.get_std('DT50_log')  # standard deviation hl
        self.full_data['DT50_log_spread'] = self.get_hl_spread()

        bmean, bstd, bmeanstd = self.get_bayesian_stats()  # bayesian stats considering LOQs
        self.full_data['DT50_log_bayesian_mean'] = bmean
        self.full_data['DT50_log_bayesian_std'] = bstd
        self.full_data['DT50_log_bayesian_mean_std'] = bmeanstd

        self.full_data['acidity_std'] = self.get_std('acidity')
        self.full_data['CEC_log_std'] = self.get_std('CEC_log')
        self.full_data['OC_log_std'] = self.get_std('OC_log')
        self.full_data['biomass_log_std'] = self.get_std('biomass_log')
        self.full_data['temperature_std'] = self.get_std('temperature')

        self.full_data.to_csv(self.full_data_tsv, sep='\t', index=False)

    def curate_halflife_data(self, from_csv=False):
        print("\n############# Half-life information curation ############# ")
        if from_csv:
            self.full_data = pd.read_csv(self.full_data_tsv, sep='\t')
            print("Existing file loaded from {}".format(self.full_data_tsv))
            return
        halflife_is_valid = []  # bool
        matching_spike = []  # bool
        halflife_model_category = []
        # lookup instead of recalculate smiles
        for index, row in self.full_data.iterrows():
            matching_spike.append(self.check_spike_consistency(row[self.smiles_name], row['spike_compound']))
            halflife_is_valid.append(self.check_comment_for_validity(row['halflife_comment']))
            model = row['halflife_model']
            if type(model) == float:
                model = ''
            if model == '':
                print("checking comment for:", row['halflife_comment'])
                model = self.check_comment_for_model_info(row['halflife_comment'])
            halflife_model_category.append(self.unify_model_descriptions(model))
        self.full_data['halflife_is_valid'] = halflife_is_valid
        self.full_data['matching_spike'] = matching_spike
        self.full_data['halflife_model_category'] = halflife_model_category

        self.full_data.to_csv(self.full_data_tsv, sep='\t', index=False)
        print('Curated file saved to', self.full_data_tsv)

    def reduce_for_modelling(self, from_csv = False): #todo: separate function to only load model data
        print("\n############# Reduce data set ############# ")
        if from_csv:
            self.cpd_data = pd.read_csv(self.cpd_data_tsv, sep='\t')
            self.model_data = pd.read_csv(self.model_data_tsv, sep='\t')
            print("Existing files loaded from {} and {}".format(self.cpd_data_tsv,self.model_data_tsv))
            return

        self.reduce_data()
        # curate manually and save again
        self.curate_manually()
        # create modelling input
        self.create_modelling_input()

    def reduce_data(self):
        # reduce dataset
        print('Data frame size: ', len(self.full_data))
        self.cpd_data = self.full_data.loc[:,
              [self.id_name, self.smiles_name, self.compound_name, 'compound_id', # 'node_depth', removed, not available anymore
               'DT50_count', 'DT50_gmean', 'DT50_log_median', 'DT50_log_gmean',
               'DT50_log_spread', 'DT50_log_std', 'DT50_log_bayesian_mean', 'DT50_log_bayesian_std',
               'DT50_log_bayesian_mean_std', 'acidity_std', 'CEC_log_std', 'OC_log_std', 'biomass_log_std', 'temperature_std',
               'canonical_SMILES', 'cropped_canonical_SMILES', 'cropped_canonical_SMILES_no_stereo']]
        self.cpd_data = self.cpd_data.drop_duplicates(self.id_name)
        self.cpd_data[self.target_variable_name] = self.cpd_data['DT50_log_bayesian_mean']
        self.cpd_data[self.target_variable_std_name] = self.cpd_data['DT50_log_bayesian_std']

        # save and describe
        print('Data frame size: ', len(self.cpd_data))
        self.cpd_data.to_csv(self.cpd_data_tsv, sep='\t', index=False)
        print('Compound data file saved to', self.cpd_data_tsv)
        description = self.cpd_data.describe()
        description.to_csv(self.cpd_data_description_file, sep='\t', index=False)

    def curate_manually(self):
        print("\n############# Manual curation ############# ")
        # Load data
        print('Original number of compounds: {}'.format(self.cpd_data.shape[0]))

        # 1. only keep organic substances
        inorganic_substances = ['Disodium phosphonate', 'Phosphine', 'Phosphorous acid', 'Zinc phosphide', 'CO2']
        for i in inorganic_substances:
            self.cpd_data = self.cpd_data.loc[self.cpd_data[self.compound_name] != i]
        print('Organic substances: {}'.format(self.cpd_data.shape[0]))

        # 2. Fix Fluroxypyr ester + acid according to information in DAR
        # change name to Fluroxypyr acid
        fluroxypyr_acid_index = int(self.cpd_data.loc[self.cpd_data[self.compound_name] == 'Fluroxypyr ester + acid'].index[0])
        self.cpd_data.at[fluroxypyr_acid_index, self.compound_name] = 'Fluroxypyr acid'
        # For Dimethomorph, there are 2 isomers and 1 mixture. Change name.
        dimethomorph_index = int(self.cpd_data.loc[self.cpd_data['compound_id'] == 'https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/e98cd963-288f-46c2-b486-6075ccd5d4ef/structure/6a08f6fd-a4fb-4568-959b-d2ada916ac9c'].index[0])
        self.cpd_data.at[dimethomorph_index, self.compound_name] = 'Dimethomorph'

        # 3. remove composite mixtures where no main substance can be determined,
        # or separate data for each compound is available
        self.cpd_data.drop(self.cpd_data[self.cpd_data[self.compound_name] == 'Guazatine'].index, inplace=True)
        self.cpd_data.drop(self.cpd_data[self.cpd_data[self.compound_name] == 'Sodium 5-nitrocompounds'].index, inplace=True)
        self.cpd_data.drop(self.cpd_data[self.cpd_data[self.compound_name] == 'combined Propaquizafop & Propaquizafop acid'].index,
                           inplace=True)
        self.cpd_data.drop(self.cpd_data[self.cpd_data[self.compound_name] == 'sum of Kresoxim-methyl and BF 490-1'].index,
                           inplace=True)

        print('Data set without composite mixtures: {}'.format(self.cpd_data.shape[0]))

        # Check:
        duplicates = self.cpd_data[self.cpd_data.duplicated([self.smiles_name])]
        print('Number of duplicated cropped canonical SMILES without stereo information:', duplicates.shape[0])
        self.cpd_data.to_csv(self.cpd_data_tsv, sep='\t', index=False)

    def check_spike_consistency(self, compound, spike):
        if type(spike) == float:
            return False

        if self.spike_compound_dictionary.get(spike):
            spike_set = self.spike_compound_dictionary[spike]
        else:
            spike_set = set()
            for species in spike.split('.'):
                clean_spike_smiles = Util.remove_isotope_info(Util.remove_stereo_info(species))
                can_spike_smiles = Util.canonicalize_smiles(clean_spike_smiles)
                spike_set.add(can_spike_smiles)
            self.spike_compound_dictionary[spike] = spike_set
        if compound in spike_set:
            return True
        else:
            return False

    def get_main_soil_texture(self):
        new = []
        for index, row in self.full_data.iterrows():
            sand = float(row['sand'])
            silt = float(row['silt'])
            clay = float(row['clay'])
            if np.isnan(sand) or np.isnan(silt) or np.isnan(clay):
                new.append('N/A')
            else:
                if sand > silt and sand > clay:
                    new.append('sand')
                elif silt > sand and silt > clay:
                    new.append('silt')
                else:
                    new.append('clay')
        return new

    def get_hl_median(self):
        new = []
        for index, row in self.full_data.iterrows():
            this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            median = np.median(this['DT50_log'])
            new.append(median)
        return new

    @staticmethod
    def g_mean(x):
        a = np.log(x)
        return np.exp(a.mean())

    def get_geometric_mean(self):
        new = []
        for index, row in self.full_data.iterrows():
            this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            gmean = self.g_mean(this['reported_DT50'])
            new.append(gmean)
        return new

    def get_std(self, column):
        new = []
        for index, row in self.full_data.iterrows():
            this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            std = np.nanstd(this[column])
            new.append(std)
        return new

    def get_hl_spread(self):
        new = []
        for index, row in self.full_data.iterrows():
            this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            spread = max(this['DT50_log']) - min(this['DT50_log'])
            new.append(spread)
        return new

    def index_compounds(self):
        new = []
        this_id = 0
        D = {}
        for index, row in self.full_data.iterrows():
            if row[self.smiles_name] not in D.keys():
                this_id += 1
                D[row[self.smiles_name]] = this_id
                new.append(this_id)
            else:
                new.append(D[row[self.smiles_name]])
        return new

    def count_halflives(self):
        new = []
        for i in self.full_data[self.id_name]:
            new.append(self.full_data[self.id_name].value_counts()[i])
        return new

    def get_bayesian_stats(self, curate_data=False):
        mean_list = []
        std_list = []
        mean_std_list = []
        results = {}  # {'index': (mean, std, mean_std)}
        for index, row in self.full_data.iterrows():
            if row[self.id_name] in results.keys():
                mean, std, mean_std = results[row[self.id_name]]
            else:
                this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
                comment_list_raw = self.process_comment_list(this["halflife_comment"])
                y_raw = np.array(this['DT50_log'])
                if curate_data:
                    is_valid = this['halflife_is_valid']
                    models = this['halflife_model_category']
                    is_spike = this['matching_spike']
                    y, comment_list = self.curate_data_points(y_raw, comment_list_raw, is_valid, models, is_spike)
                else:
                    y = y_raw
                    comment_list = comment_list_raw
                print("\nCOMPOUND INDEX {}".format(row[self.id_name]))
                print("Compute bayes for {} with comments {}".format(y, comment_list))
                bayesian = Bayesian(y=y, comment_list=comment_list)
                bayesian.set_prior_mu(mean=1.5, std=2)
                bayesian.set_prior_sigma(mean=0.4, std=0.4)
                bayesian.set_lower_limit_sigma(0.2)
                mean, std, mean_std = bayesian.get_posterior_distribution()
                results[row[self.id_name]] = (mean, std, mean_std)
                print('mean: {}, std: {}, mean_std: {}'.format(mean, std, mean_std))
            mean_list.append(round(mean, 2))
            std_list.append(round(std, 2))
            mean_std_list.append(round(mean_std, 2))
        return mean_list, std_list, mean_std_list

    @staticmethod
    def process_comment_list(comment_list):
        new_list = []
        for comment in comment_list:
            if type(comment) == float:
                new_list.append('')
            elif '<' in comment:
                new_list.append('<')
            elif '>' in comment:
                new_list.append('>')
            else:
                new_list.append('')
        return new_list

    def curate_data_points(self, y, comment, is_valid, models, is_spike):
        remove_indexes = []
        # remove invalid
        for index, v in enumerate(is_valid):
            if v == False and True in is_valid:
                remove_indexes.append(index)
        y, comment, is_valid, models, is_spike = self.remove_indexes_from_list(y, comment, is_valid, models, is_spike,
                                                                               remove_indexes)

        # remove non-SFO
        remove_indexes = []
        for index, m in enumerate(models):
            if 'SFO' in models and 'First Order' in models:
                if m not in ['SFO', 'First Order']:
                    remove_indexes.append(index)
        y, comment, is_valid, models, is_spike = self.remove_indexes_from_list(y, comment, is_valid, models, is_spike,
                                                                               remove_indexes)

        # remove non matching spike
        remove_indexes = []
        for index, s in enumerate(is_spike):
            if s == True and False in is_spike:
                remove_indexes.append(index)
        y, comment, is_valid, models, is_spike = self.remove_indexes_from_list(y, comment, is_valid, models, is_spike,
                                                                               remove_indexes)
        return y, comment

    @staticmethod
    def remove_indexes_from_list(list1, list2, list3, list4, list5, indexes):  # ok this is veery ugly
        output_list1 = []
        output_list2 = []
        output_list3 = []
        output_list4 = []
        output_list5 = []
        for index, item in enumerate(list1):
            if index not in indexes:
                output_list1.append(item)
        for index, item in enumerate(list2):
            if index not in indexes:
                output_list2.append(item)
        for index, item in enumerate(list3):
            if index not in indexes:
                output_list3.append(item)
        for index, item in enumerate(list4):
            if index not in indexes:
                output_list4.append(item)
        for index, item in enumerate(list5):
            if index not in indexes:
                output_list5.append(item)
        return output_list1, output_list2, output_list3, output_list4, output_list5

    @staticmethod
    def check_comment_for_validity(comment):
        if type(comment) == float:
            comment = ''
        valid = True
        for search_tag in ["not acceptable", "not acceptabl", "not considered valid",
                           "should not be used as imput for modelling",
                           "not considered representative for environmental conditions", "implausible",
                           "RMS recommends not using data from this study for fate and behaviour assessment"]:
            if search_tag in comment:
                valid = False
        return valid

    @staticmethod
    def check_comment_for_model_info(comment):
        if type(comment) == float:
            comment = ''
        comment_lower = comment.lower()
        model = ''
        for search_tag in ["SFO", "FOMC", "DFOP", "first-order", "first order", "FOMC-DFOP", "FMOC", "HS",
                           "Hockey stick", "hockey-stick", "SFO-SFO"]:
            if search_tag.lower() in comment_lower:
                model = search_tag
        return model

    @staticmethod
    def unify_model_descriptions(model):
        new = str(model)
        if model == "FMOC":
            new = "FOMC"
        elif model.lower() in ["hockey stick", "hockey-stick"]:
            new = "HS"
        elif model.lower() in ["first-order", "first order", "first order kinetics"]:
            new = "First Order"
        elif model == "FOMC + SFO":
            new = "FOMC-SFO"
        elif model in ["SFO/SFO", "SFO + SFO"]:
            new = "SFO-SFO"
        elif model == "SFOÂ²":
            new = "SFO"
        return new

    # Visualization
    def analyze_distributions(self):
        """
        This function visualizes the distribution of the mean and the standard deviation of the target variable.
        Optionally, distributions obtained from Bayesian inference can be considered.
        """
        print("\n############# Analyze target variable distribution ############# ")
        v = Visualize(self,'analyze_distributions')
        # distribution of target variable
        v.plot_target_variable_distribution(mean_name='DT50_log_gmean', std_name='DT50_log_std', cutoff_value = 20,
                                            include_BI = True,
                                          BI_mean_name = 'DT50_log_bayesian_mean',
                                          BI_std_name = 'DT50_log_bayesian_std')

