from pepper_lab.pepper import Pepper
from pepper_lab.util import *
from pepper_lab.datastructure import DataStructure
from pepper_lab.bayesian import *


class DataStructureSediment(DataStructure):
    def __init__(self, pep: Pepper):
        super().__init__(pep)
        self.set_data_directory(os.data.join(pep.data_directory,'data_structure','sediment'))
        self.envipath_package = 'https://envipath.org/package/5c5639b0-19f8-4671-9d9a-36f05e5518df'  # bigger
        # self.envipath_package = 'https://envipath.org/package/833d620c-db01-4650-9ef8-2c5de2edf3dd'  # smaller
        self.data_type = 'sediment'
        self.data_dict = {}
        self.raw_data_tsv = self.build_output_filename('raw_data')
        self.full_data_tsv = self.build_output_filename('full_data')
        self.cpd_data_tsv = self.build_output_filename('cpd_data')
        self.model_data_tsv = self.build_output_filename('model_data')
        self.cpd_data_description_file = self.build_output_filename('cpd_data_description')
        self.spike_compound_dictionary = {}

    def curate_annotate(self, from_csv: bool = False, from_paper: bool = False):
        """
        Curate SMILES, log-transform values, calculate bayesian-inferred mean of target variable,
        and curate information on half-lives etc.
        :param from_csv: If true, load existing csv files
        :param from_paper: If true, load existing csv file
        :return:
        """
        if from_csv:
            if from_paper:
                raise NotImplementedError
            else:
                self.full_data = pd.read_csv(self.full_data_tsv, sep='\t')
                print("Existing curated file loaded from {}".format(self.full_data_tsv))
                return
        else:
            self.full_data = self.raw_data

        self.curate_smiles()
        self.transform_values()
        self.curate_halflife_data()

    def transform_values(self, from_csv=False):
        """
        This function extracts sediment texture,study_guideline information
        and log-transform values, computes the bayesian-inferred mean of target variable etc.
        :param from_csv: If true, load existing csv file
        :return:
        """
        print("\n############# Value transformation - log transformation and bayesian inference ############# ")
        if from_csv:
            self.full_data = pd.read_csv(self.full_data_tsv, sep='\t')
            print("Existing file loaded from {}".format(self.full_data_tsv))
            return

        # Index compounds by SMILES identity
        self.full_data[self.id_name] = self.index_compounds()
        self.full_data.sort_values(by=self.id_name)

        # Log values calculation
        self.full_data['DT50_log_total_system'] = Util.log_transform(self.full_data['DT50_total_system'])
        self.full_data['DT50_log_water'] = Util.log_transform(self.full_data['DT50_water'])
        self.full_data['DT50_log_sediment'] = Util.log_transform(self.full_data['DT50_sediment'])

        self.full_data['soil_texture_main'] = self.get_main_soil_texture()
        self.full_data['sediment_texture_main'] = self.get_main_sediment_texture()
        self.full_data['study_guideline'] = self.get_study_guideline()

        self.full_data['CEC_log'] = Util.log_transform(self.full_data['CEC'])
        self.full_data['biomass_log'] = Util.log_transform(self.full_data['biomass'])
        self.full_data['OC_log'] = Util.log_transform(self.full_data['OC'])

        self.full_data['OM_log'] = Util.log_transform(self.full_data['OM'])
        self.full_data['TOC_log'] = Util.log_transform(self.full_data['TOC'])
        self.full_data['DOC_log'] = Util.log_transform(self.full_data['DOC'])

        self.full_data['DT50_count'] = self.count_halflives_ws()

        self.full_data['DT50_total_system_log_mean'] = self.get_hl_mean_ws()  # half-lives (hl) log mean
        self.full_data['DT50_total_system_log_median'] = self.get_hl_median_ws()  # log median
        self.full_data['DT50_log_std_total_system'] = self.get_std('DT50_log_total_system')  # standard deviation hl
        self.full_data['DT50_log_spread'] = self.get_hl_spread_ws()

        b_mean, b_std, b_mean_std = self.get_bayesian_stats()  # bayesian stats
        self.full_data['DT50_log_bayesian_mean'] = b_mean
        self.full_data['DT50_log_bayesian_std'] = b_std
        self.full_data['DT50_log_bayesian_mean_std'] = b_mean_std

        self.full_data.to_csv(self.full_data_tsv, sep='\t', index=True)

    def curate_halflife_data(self, from_csv=False):
        """
        This function curates the information on half-lives.
        :param from_csv: If true, load existing csv file
        :return:
        """
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

    def reduce_for_modelling(self, from_csv=False):
        print("\n############# Reduce data set ############# ")
        if from_csv:
            self.cpd_data = pd.read_csv(self.cpd_data_tsv, sep='\t')
            self.model_data = pd.read_csv(self.model_data_tsv, sep='\t')
            print("Existing files loaded from {} and {}".format(self.cpd_data_tsv, self.model_data_tsv))
            return

        self.reduce_data()
        # create modelling input
        self.create_modelling_input()

    def reduce_data(self):
        """
        Compiles the relevant data properties for each compound and its corresponding data summary
        of all the properties.
        :return: saves two csv files, namely, compound data file and compound description (data summary) file
        """
        print('Data frame size: ', len(self.full_data))
        self.cpd_data = self.full_data.loc[:,
                        [self.id_name, self.smiles_name, 'compound_name', 'compound_id', 'DT50_count',
                         'DT50_total_system_log_mean', 'DT50_total_system_log_median',
                         'DT50_log_spread', 'DT50_log_std_total_system', 'DT50_log_bayesian_mean',
                         'DT50_log_bayesian_std', 'DT50_log_bayesian_mean_std',
                         'canonical_SMILES', 'cropped_canonical_SMILES', 'cropped_canonical_SMILES_no_stereo']]
        self.cpd_data = self.cpd_data.drop_duplicates(self.id_name)
        self.cpd_data[self.target_variable_name] = self.cpd_data['DT50_log_bayesian_mean']
        self.cpd_data[self.target_variable_std_name] = self.cpd_data['DT50_log_bayesian_std']

        # save compound data and data summary
        print('Data frame size: ', len(self.cpd_data))
        self.cpd_data.to_csv(self.cpd_data_tsv, sep='\t', index=False)  # compounds' data
        data_description = self.cpd_data.describe()  # data summary of key properties
        data_description.to_csv(self.cpd_data_description_file, sep='\t', index=False)

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

    def get_main_sediment_texture(self):
        new = []
        for index, row in self.full_data.iterrows():
            sand = float(row['sand'])
            silt = float(row['silt'])
            clay = float(row['clay'])
            if np.isnan(sand) or np.isnan(silt) or np.isnan(clay):
                new.append('N/A')
            else:
                if clay + silt > 50:
                    new.append('Fine texture')
                elif clay + silt < 50:
                    new.append('Coarse texture')
                elif clay + silt == 50:
                    new.append('borderline')
                else:
                    new.append('')
        return new

    def get_study_guideline(self):
        study_name = []
        for index, row in self.full_data.iterrows():
            study_description = row['study_description']
            if 'OECD308' in study_description or 'OECD 308' in study_description:
                study_name.append('OECD 308')
            elif 'SETAC' in study_description or 'Setac' in study_description:
                study_name.append('SETAC')
            elif 'BBA' in study_description:
                study_name.append('BBA')
            elif 'EPA' in study_description:
                study_name.append('EPA')
            else:
                study_name.append('others')

        return study_name

    def get_hl_mean_ws(self):
        new = []
        for index, row in self.full_data.iterrows():
            this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            mean = np.mean(this['DT50_log_total_system'])
            new.append(mean)
        return new

    def get_hl_median_ws(self):
        new = []
        for index, row in self.full_data.iterrows():
            this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            median = np.median(this['DT50_log_total_system'])
            new.append(median)
        return new

    def get_std(self, column):
        new = []
        for index, row in self.full_data.iterrows():
            this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            std = np.nanstd(this[column])
            new.append(std)
        return new

    def get_hl_spread_ws(self):
        new = []
        for index, row in self.full_data.iterrows():
            this = self.full_data.loc[self.full_data[self.id_name] == row[self.id_name]]
            spread = max(this['DT50_log_total_system']) - min(this['DT50_log_total_system'])
            new.append(spread)
        return new

    def index_compounds(self):
        new = []
        this_id = 0
        d = {}
        for index, row in self.full_data.iterrows():
            if row[self.smiles_name] not in d.keys():
                this_id += 1
                d[row[self.smiles_name]] = this_id
            new.append(this_id)
        return new

    def count_halflives_ws(self):
        new = []
        for i in self.full_data[self.id_name]:
            new.append(self.full_data[self.id_name].value_counts()[i])
        return new

    def get_bayesian_stats(self):
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
                y_raw = np.array(this['DT50_log_total_system'])
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

    # todo: Visualization
    # def visualize(self):
    #     vis = VisualizeSediment(self.pepper, self)
    #     vis.print_violin_plot()

    @staticmethod
    def check_for_kinetics(addinfo):
        try:
            addinfo.get_halflife_ws().get_value()
        except AttributeError:
            return False
        else:
            return True
