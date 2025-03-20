from tqdm import tqdm

from pepper.pepper import Pepper
from pepper.util import *
from pepper.datastructure import DataStructure
from pepper.bayesian import *


class DataStructureSludge(DataStructure):
    def __init__(self, pep: Pepper):
        super().__init__(pep)
        self.set_data_directory(os.path.join(pep.data_directory,'data_structure','sludge'))
        self.smiles_name = pep.get_smiles_name()
        self.target_variable_name = pep.get_target_variable_name()
        self.target_variable_std_name = pep.get_target_variable_std_name()
        self.compound_name = pep.get_compound_name()
        self.id_name = pep.get_id_name()

        self.envipath_package = 'https://envipath.org/package/7932e576-03c7-4106-819d-fe80dc605b8a'
        self.data_type = 'sludge'
        self.data_dict = {}

        # These are general, yet, they must be initialized with the info of the new object
        # Can we avoid defining all these for every subclass?
        self.raw_data_tsv = self.build_output_filename('raw_data')
        self.full_data_tsv = self.build_output_filename('full_data')
        self.cpd_data_tsv = self.build_output_filename('cpd_data')
        self.model_data_tsv = self.build_output_filename('model_data')
        self.cpd_data_description_file = self.build_output_filename('cpd_data_description')

    def curate_annotate(self, from_csv: bool = False):
        """
        have either rate constants or half-lives, we can utilize the reaction order formula to convert kinetic data
        between them. The default reaction type assumed is a first-order reaction. In addition to this, the script also
        incorporates biomass for each half-life. Finally, the relevant kinetic data is collected and subjected to log
        calculations, including

        :param from_csv: If true, load existing csv files
        :return:
        """
        if from_csv:
            self.full_data = pd.read_csv(os.path.join('..', 'data', 'full_data_sludge_all_data.tsv'), sep='\t')
            return

        self.full_data = self.raw_data
        self.curate_smiles()
        self.transform_values()

    def transform_values(self):
        self.full_data["kinetics_comment"] = self.full_data[["rateconstant_comment", "halflife_comment"]].apply(lambda x: self.process_comment_list(x), axis=1)
        self.full_data['k_combined'] = self.full_data.apply(lambda x: self.get_k(x), axis=1)
        self.full_data['k_biomass_corrected'] = self.full_data.apply(lambda x: self.get_k_biomass(x), axis=1)
        self.full_data['halflife'] = self.full_data.apply(lambda x: self.get_DT50(x), axis=1)
        self.full_data['hl_biomass_corrected'] = self.full_data.apply(lambda x: self.get_DT50_biomass(x), axis=1)
        self.full_data['log_k_combined'] = np.log10(self.full_data['k_combined'])
        self.full_data['log_k_biomass_corrected'] = np.log10(self.full_data['k_biomass_corrected'])
        self.full_data['halflife_log'] = np.log10(self.full_data['halflife'])
        self.full_data['log_hl_biomass_corrected'] = np.log10(self.full_data['hl_biomass_corrected'])

        self.calculate_target_variables()
        self.full_data.dropna(subset=['halflife', 'halflife_log'], inplace=True)
        self.calculate_bay_mean_std()
        self.full_data.to_csv(self.full_data_tsv, sep='\t', index=False)


    def reduce_for_modelling(self, from_csv = False):
        print("\n############# Reduce data set ############# ")
        if from_csv:
            self.cpd_data = pd.read_csv(self.cpd_data_tsv, sep='\t')
            self.model_data = pd.read_csv(self.model_data_tsv, sep='\t')
            print("Existing files loaded from {} and {}".format(self.cpd_data_tsv,self.model_data_tsv))
            return

        self.reduce_data()
        # create modelling input
        self.create_modelling_input()

    def calculate_target_variables(self):
        self.full_data[['hl_gmean', 'biomass_hl_gmean']] = self.full_data.groupby('canonical_SMILES')[['halflife', 'hl_biomass_corrected']].transform(lambda x: self.g_mean(x))
        self.full_data[['hl_log_median', 'biomass_hl_log_median']] = self.full_data.groupby('canonical_SMILES')[['halflife_log', 'log_hl_biomass_corrected']].transform('median')
        self.full_data[['hl_log_std', 'biomass_hl_log_std']] = self.full_data.groupby('canonical_SMILES')[['halflife_log', 'log_hl_biomass_corrected']].transform(lambda x: np.nanstd(x))
        self.full_data[['acidity_std', 'temperature_std',
                        'biomass_std']] = self.full_data.groupby('canonical_SMILES')[['acidity', 'temperature', 'total_suspended_solids_concentration_start']].transform(lambda x: np.nanstd(x))
        self.full_data[['hl_log_spread', 'biomass_hl_log_spread']] = self.full_data.groupby('canonical_SMILES')[['halflife_log', 'log_hl_biomass_corrected']].transform(lambda x: max(x) - min(x))
        return self.full_data

    def reduce_data(self):
        # reduce dataset
        print('Data frame size: ', len(self.full_data))
        self.cpd_data = self.full_data.loc[:, [
            self.id_name, self.smiles_name, self.compound_name, 'compound_id',
            'hl_gmean', 'biomass_hl_gmean', 'hl_log_median', 'biomass_hl_log_median',
            'hl_log_std', 'biomass_hl_log_std', 'hl_log_spread', 'biomass_hl_log_spread',
            'hl_log_bayesian_mean', 'hl_log_bayesian_std',
            'biomass_std', 'acidity_std', 'temperature_std',
            'canonical_SMILES', 'cropped_canonical_SMILES', 'cropped_canonical_SMILES_no_stereo']]
        self.cpd_data = self.cpd_data.drop_duplicates(self.id_name)
        self.cpd_data.rename(columns={"DT50_log_bayesian_mean": self.target_variable_name,
                                      "DT50_log_bayesian_std": self.target_variable_std_name},
                             inplace=True)

        # save and describe
        print('Data frame size: ', len(self.cpd_data))
        self.cpd_data.to_csv(self.cpd_data_tsv, sep='\t', index=False)
        description = self.cpd_data.describe()
        description.to_csv(self.cpd_data_description_file, sep='\t', index=False)

    def process_comment_list(self, df):
        if df.str.contains("<").any() or df.str.contains(">").any():
            # This line adjusts the fact that rateconstant and halflifes are inversely related
            # i.e. a lower truncated rateconstant represents an upper truncated halflife
            if pd.Series([df.str.contains(">")[0], df.str.contains("<")[1]]).any():
                return '<'
            else:
                return '>'
        else:
            return ''

    def get_k(self, row):
        k_given = row['rateconstant']
        k_unit = row['rateconstant_unit']
        k_true = np.NaN
        TSS = row['total_suspended_solids_concentration_start']
        hl = row['halflife_raw']
        order = row['halflife_model']
        if not np.isnan(k_given):
            if k_given != 0 and k_unit == '1 / day' and not np.isnan(TSS):
                k_true = k_given
            elif k_given != 0 and k_unit == 'L / (g TSS * day)' and not np.isnan(TSS):
                k_true = k_given * TSS
            elif k_given != 0 and k_unit == 'μg / (g TSS * day)' and not np.isnan(TSS): #TODO: Should this be merged with above?
                pass
            else:
                if np.isnan(TSS):
                    print(f'Problem: no TSS for scenario {row["scenario_id"]}')
                elif k_given == 0:
                    print(f'Problem: given rate constant is 0 for scenario {row["scenario_id"]}')
        elif not np.isnan(hl):
            if order == 'Zero order':
                k_true = TSS / (2 * hl)
            elif order == 'First order':
                k_true = np.log(2) / hl
            elif order == 'Pseudo first order':  # it's a biomass corrected hl
                real_hl = hl / TSS
                k_true = np.log(2) / real_hl
            else:
                k_true = np.log(2) / hl  # By default, using the 1st order reaction formula
        else:
            print(f'rate constant is na and halflife too for scenario {row["scenario_id"]}')
        return k_true

    def get_k_biomass(self, row):
        k_given = row['rateconstant']
        k_unit = row['rateconstant_unit']
        TSS = row['total_suspended_solids_concentration_start']
        hl = row['halflife_raw']
        k_biomass = np.NaN
        k = row['k_combined']
        if not np.isnan(k_given) and k_given != 0:
            if k_unit == '1 / day':
                k_biomass = k_given / TSS
            elif k_unit == 'L / (g TSS * day)': #TODO: Add 'μg / (g TSS * day)'?
                k_biomass = k_given
        elif k_given == 0:
            print(f'Error: rate constant is 0 for scenario {row["scenario_id"]}')
        elif np.isnan(k_given) and not np.isnan(hl):
            k_biomass = k / TSS
        return k_biomass

    def get_DT50(self, row):
        k = row['k_combined']
        hl_given = row['halflife_raw']
        hl = np.NaN
        if np.isnan(hl_given) and not np.isnan(k):
            hl = np.log(2) / k
        elif hl_given != 0:
            hl = hl_given
        else:
            print(f'Error: half-life == 0 for scenario {row["scenario_id"]}')
        return hl

    def get_DT50_biomass(self, row):
        TSS = row['total_suspended_solids_concentration_start']
        hl = row['halflife']
        k_biomass = row['k_biomass_corrected']
        hl_biomass = np.NaN
        if not np.isnan(hl):
            hl_biomass = hl / TSS
        elif not np.isnan(k_biomass):
            hl_biomass = np.log(2) / k_biomass
        return hl_biomass

    def g_mean(self, x):
        a = np.log(x)
        return np.exp(a.mean())

    def calculate_bay_mean_std(self):
        bmean, bstd = self.get_bayesian_stats()
        self.full_data['hl_log_bayesian_mean'] = bmean
        self.full_data['hl_log_bayesian_std'] = bstd
        # df.to_csv(output_file_path+'sludge_calculated_test_for_baycalculation.tsv', sep='\t')
        return

    def get_bayesian_stats(self):
        mean_list = []
        std_list = []
        results = {}  # {'index': (mean, std)}
        for index, row in tqdm(self.full_data.iterrows()):
            if row['canonical_SMILES'] in results.keys():
                mean, std = results[row['canonical_SMILES']]
            else:
                df = self.full_data.loc[self.full_data['canonical_SMILES'] == row['canonical_SMILES']]
                y = np.array(df['halflife_log'])
                kinetics_comments = df.kinetics_comment.to_list()
                print("\nCOMPOUND canonical_SMILES {}".format(row['canonical_SMILES']))
                print("Compute bayes for {} with comments {}".format(y, kinetics_comments))
                bayesian = Bayesian(y=y, comment_list=kinetics_comments)
                bayesian.set_prior_mu(mean=0, std=1)  # (Original: mean=1.5, std=2) Set prior_mu_std as 2
                bayesian.set_prior_sigma(mean=0.6, std=0.2)  # (Original: mean=0.2, std=0.5)
                bayesian.set_lower_limit_sigma(0.3)
                mean, std, _ = bayesian.get_posterior_distribution()
                results[row['canonical_SMILES']] = (mean, std)
                print('mean: {}, std: {}'.format(mean, std))
            mean_list.append(round(mean, 2))
            std_list.append(round(std, 2))

        return mean_list, std_list
