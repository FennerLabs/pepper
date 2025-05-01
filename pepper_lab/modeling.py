from fileinput import filename

import pandas as pd
import numpy as np
import os
import sys
from copy import deepcopy
import yaml

from pepper_lab.datastructure import DataStructure
from pepper_lab.pepper import Pepper
from pepper_lab.descriptors import Descriptors
from pepper_lab.model import Model
from pepper_lab.visualize import Visualize
from pepper_lab.util import Util


# exploratory analysis, training, evaluation
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline


# Importing all regressors
# from sklearn.linear_model import LinearRegression
# from sklearn.linear_model import Ridge
# from sklearn.linear_model import SGDRegressor
# from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.neural_network import MLPRegressor
# from sklearn.tree import DecisionTreeRegressor


class Modeling(Pepper):
    def __init__(self, pep: Pepper(), data: DataStructure, descriptors: Descriptors):
        super().__init__()
        self.set_data_directory(os.path.join(pep.data_directory, 'modeling'))
        self.tag = pep.get_tag()
        self.data_type = pep.get_data_type()
        self.curation_type = pep.get_curation_type()
        self.target_variable_name = pep.get_target_variable_name()
        self.target_variable_std_name = pep.get_target_variable_std_name()
        self.compound_name = pep.get_compound_name()
        self.id_name = pep.get_id_name()
        self.smiles_name = pep.get_smiles_name()
        self.pepper = pep
        self.random_state = pep.get_random_state()

        self.reduced_features = pd.DataFrame()  # subset of features after feature reduction

        # Dataframe including both features and endpoint
        self.joint_data = pd.DataFrame()

        # Attributes from data that we want to keep
        self.model_data = data.model_data
        self.descriptors = descriptors
        self.target_variable = data.get_target_variable()

        # Defining the default pipe
        self.variance_selector = VarianceThreshold()
        self.scaler = MinMaxScaler()
        self.regressor = RandomForestRegressor(random_state=0, n_jobs=-1)

        # Scoring storage
        self.test_scores = pd.DataFrame()
        self.train_scores = pd.DataFrame()
        self.train_scores_tsv = None
        self.test_scores_tsv = None
        self.predicted_target_variable_tsv = None

        # feature space to be explored
        self.feature_space_list = [] # e.g., ['maccs', 'maccs+padel', 'all']
        self.complete_feature_space_list = descriptors.get_feature_space_list()

        # Regressor settings
        # all available regressor settings imported from yaml config file: {'RF': {'regressor':RandomForest(), ...}}
        self.regressor_settings = self.load_regressor_settings()
        # This is a customized list of dictionaries from regressor_settings that are used for modelling
        self.regressor_list = [] # [{'regressor': RandomForest(), ...}]
        # List of regressors names used
        self.regressor_name_list = [] #e.g., ['RF', 'SVR', ...]
        # List of all regressor names available within PEPPER
        self.complete_regressor_name_list = []

        # Best models from
        self.best_models = []



    #-----------------------------------------------------#
    # Key functions, includes those with cross validation #
    #-----------------------------------------------------#

    def run_test_model(self, regressor_name = 'RF'):
        """
        # Run a test model to verify correct data and descriptor config
        """
        print("############# Run test model #############")
        # get the name of the function
        function_name = sys._getframe().f_code.co_name
        # define which descriptors to use
        self.descriptors.define_feature_space('all')
        # define model object
        test_model = Model(self.pepper, self.descriptors, self.model_data)
        # prepare X and y for modeling
        test_model.prepare()
        # train model on a single split and evaluate
        test_model.simple_run()
        # collect scores from simple_run
        self.collect_scores(test_model)
        # visualize results
        v = Visualize(test_model, function_name)
        v.scatterplot_predicted_vs_test()

    def test_different_setups_cross_val(self, number_of_splits=5,
                                        regressor_name_list=None, feature_space_list=None,
                                        include_plant_fingerprints=False, load_existing=False, config='default', **kwargs):
        """
        Test different combinations of regressors and descriptors with k-fold cross-validation.
        The different combinations are created with the method 'setting_the_stage'.
        First a reference model is created for each (regressor, descriptor) combination.
        Then there is an outer loop based on k-fold cross validation (i.e., "outer loop").
        This is simpler than the 'nested_cross_validation()' method which later splits the data again
        in an 'inner loop' for optimization.

        :param regressor_name_list: List of regressors to be compared.
        If none set, an extensive list of regressors will be used.
        :param feature_space_list: List fo feature spaces to be compared
        :param include_plant_fingerprints: Default behaviour is to ignore this.
        Currently used to account cases in which the plant  conditions are described using a 'plant fingerprint'.
        :param number_of_splits: Number of splits for k-fold cross-validation
        :param load_existing: Load existing results
        :param config: regressor configs as defined in ./config/regressor_settings__['range'/'singlevalue']_[config].yml config file
        """

        print("\n############# Test different setups cross-validation #############")
        function_name = sys._getframe().f_code.co_name  # get the name of the function

        # The regressor_name_list is created here by calling the setting_the_stage method
        self.setting_the_stage(function_name, regressor_name_list, feature_space_list,
                               mode='singlevalue', config=config, load_existing=load_existing)

        # Iterating through feature spaces to be explored
        for feature_space in self.feature_space_list:
            self.descriptors.define_feature_space(feature_space)
            reference_model = self.prepare_reference_model()
            # Create a KFold object with n_splits
            kf = KFold(n_splits=number_of_splits, shuffle=True, random_state=self.random_state)
            # Iterate through regressor list
            for regressor_dict in self.regressor_list:
                reference_model.load_settings_from_dict(regressor_dict)
                # no feature slection/reduction is performed in this screening
                reference_model.feature_selection_method = None
                reference_model.feature_reduction_method = None
                # Define the "outer loop"
                fold_number = 1
                for train_index, test_index in kf.split(reference_model.data):

                    print(f"{'=' * 50}\n{regressor_dict['name']}, CV {fold_number}\n{'=' * 50}")
                    # a new model object is initiated every time a new test set is defined
                    model = deepcopy(reference_model)

                    if include_plant_fingerprints:
                        model.include_plant_fingerprints()

                    # Data is split according to the "outer loop" for a typical nested cross validation workflow
                    model.define_outer_loop(train_index, test_index)

                    # The model is trained and evaluated based on the train and test sets
                    # defined on the 'define_outer_loop()' method
                    model.simple_cross_val_run(run_id=fold_number, verbose=True)

                    self.collect_scores(model)
                    fold_number += 1  # Increment the fold number after each iteration but not for each regressor
                # -----------------------------------------------------------#
                # -----------------------------------------------------------#

        self.visualize_test_scores(analysis_type=function_name, categories=['descriptors','regressor'])  #
        # print test score overview
        self.save_test_score_overview()


    def nested_cross_val_screening(self, regressor_name_list = ['RF', 'GB','AB', 'SVR', 'GPR', 'KNN'],
                                   feature_space_list = ['all'],
                                   load_existing=False, config='default'):
        """
        This method intends to easily evaluate many combinations of regressors and descriptors
        to understand which (features,regressor) pairs perform better for a given dataset.
        It is done in a typical nested cross validation fashion.
        Both "inner" and "outer" loops are 5-fold-CV

        :param regressor_name_list: List of abbreviations of regressors for which nested CV will be performed
        :param feature_space_list: Features to be considered. By default, are loaded features are considered
        :param load_existing: Load existing score files or visualisation
        :param config: regressor configs as defined in ./config/regressor_settings__['range'/'singlevalue']_[config].yml config file
        """
        function_name = sys._getframe().f_code.co_name  # get the name of the function
        print("\n############# Nested cross-validation screening #############")
        self.setting_the_stage(function_name, regressor_name_list=regressor_name_list, feature_space_list=feature_space_list,
                               mode='range', config=config, load_existing=load_existing)

        for feature_space in feature_space_list:
            self.descriptors.define_feature_space(feature_space)

            # Create reference Model object, repare data and feature matrices, preprocess features
            reference_model = self.prepare_reference_model()

            # Create a KFold object with n_splits
            kf = KFold(n_splits=5, shuffle=True, random_state=self.random_state)

            # perform nested CV with feature selection, hyperparameter tuning and evaluation on test set for each
            # combination of regressor and feature selection method
            for regressor_dict in self.regressor_list:
                # Get the settings of the regressor
                reference_model.load_settings_from_dict(regressor_dict)
                # Nested cross validation "outer loop"
                fold_number = 0
                for train_index, test_index in kf.split(reference_model.data):
                    # Increment the outer CV fold number after each iteration
                    fold_number += 1
                    print(f"{'=' * 50}\n{regressor_dict['name']}, CV {fold_number}\n{'=' * 50}")

                    # A new model object is initiated as a copy from the reference model
                    model = deepcopy(reference_model)

                    # Set config name to distinguish the different models in the output
                    model.set_setup_name('CV{}'.format(fold_number))

                    # Define the "outer loop" for a typical nested cross validation workflow
                    model.define_outer_loop(train_index, test_index)

                    # Optimize using cross validation only on the training set
                    # This is where the "inner loops" occur
                    model.complete_train_regressors()

                    # The optimized model is tested on the test set which was defined in the "outer loop"
                    model.complete_evaluate_models(run_id=fold_number)

                    # Best settings are saved and performance scores collected
                    self.best_models.append((model.regressor_name, model.best_regressor_model))
                    self.collect_scores(model)

        # Visualisation of performance scores of the optimized models on the respective test sets
        self.visualize_test_scores(analysis_type=function_name,
                                   categories=['regressor', 'descriptors', 'feature_selection', 'feature_reduction'])


    # -----------------------------------------------------#
    # Functions for optimized regressors  #
    # -----------------------------------------------------#
    def build_final_model(self, regressor_name: str, feature_space: str, config):
        """
        Build final model using optimized configurations
        :param config: file tag for config file with optimized regressor settings
        :return: trained Model() object
        """
        print("\n############# build final, optimized model using all data #############")
        function_name = sys._getframe().f_code.co_name  # get the name of the function
        self.setting_the_stage(function_name, regressor_name_list = [regressor_name],
                               feature_space_list=[feature_space], mode='singlevalue', config=config)
        regressor_dict = self.get_optimized_regressor_dict(config)
        # Create new model object
        model = Model(self.pepper, self.descriptors, self.model_data)
        # load settings from config
        model.load_settings_from_dict(regressor_dict)
        # set feature space
        model.descriptors.define_feature_space(feature_space)
        # prepare model
        model.prepare()
        # define training set
        model.use_all_data_for_training()
        # feature selection
        model.select_features()    # Here features are truly selected model.features_names -> model.selected_features_names (i.e. feature importance)
        model.reduce_feature_space() # Here the feature space is reduced selected_features_names -> reduced_features (i.e. PCA)

        # visualize selected features where applicable
        if model.feature_selection_method in ['sequential', 'importance', 'pca', 'svd']:
            v = Visualize(model, function_name)
            v.feature_selection_plot()

        # train final model()
        model.set_regressor_parameters(model.regressor, model.regressor_params)
        model.train_model()
        return model


    def test_performance_vs_size(self, regressor_name_list=None, feature_space_list=None,
                                 number_of_splits: int = 100, load_existing: bool = False):
        """
        #todo @jose: can we remove this function and just use test_performance_vs_size_optimized_setup
        To evaluate the impact of the size of the training data set, all combinations of regressors x features are used
        to train a model on 10, 20, ... to 100% of the training data.

        @param regressor_list: list of regressors to be tested. An exhaustive list of regressors used if the parameter
        is left empty. For a short list of popular regressors (Random Forest, Gradient Boost, SVR), set regressor_name_list='short_list'.
        @param feature_space_list: List of features to be used. Default
        @param number_of_splits: Number of data splits to evaluate model performance
        @param load_existing: if True, loads existing test scores for visualisation
        @return:
        """
        print("\n############# Test performance vs size #############")
        function_name = sys._getframe().f_code.co_name  # get the name of the function

        self.setting_the_stage(function_name,  regressor_name_list,  load_existing)

        # Iterate through feature spaces
        for feature_space in feature_space_list:
            self.descriptors.define_feature_space(feature_space)

            for regressor_dict in self.regressor_list:
                model = self.create_model_from_dict(regressor_dict)
                # If problems an encounter check the "reference model" solution

                # This is where the method differs from other similar methods
                # -----------------------------------------------------------#
                model.size_dependent_run(random_state_list=list(range(0, number_of_splits)),)
                # -----------------------------------------------------------#

                self.collect_scores(model)

        self.visualize_test_scores(analysis_type=function_name, categories = ['descriptors','regressor'])  #


    def test_performance_vs_size_optimized_setup(self, regressor_name: str, feature_space: str, config: str,
                                                 load_existing: bool = False):
        """
        To evaluate the impact of the size of the training data set, all combinations of regressors x features are used
        to train a model on 10, 20, ... to 100% of the training data.
        :param regressor_name: Short name of regressor to be tested
        :param feature_space: Features to be used
        :param config: dictionary defining regressor, parameters, and feature selection method
        :param load_existing: if True, loads existing test scores for visualisation
        """
        print("\n############# Test performance vs size with optimized model #############")
        function_name = sys._getframe().f_code.co_name  # get the name of the function
        self.setting_the_stage(function_name, regressor_name_list = [regressor_name],
                               feature_space_list=[feature_space], mode='singlevalue', config=config)
        regressor_dict = self.get_optimized_regressor_dict(config)
        # Create new model object
        reference_model = Model(self.pepper, self.descriptors, self.model_data)
        # load settings from config
        reference_model.load_settings_from_dict(regressor_dict)
        # set feature space
        reference_model.descriptors.define_feature_space(feature_space)
        # prepare model
        reference_model.prepare()


        # Create a KFold object with n_splits
        kf = KFold(n_splits=5, shuffle=True, random_state=self.random_state)

        fold_number = 0
        for train_index, test_index in kf.split(reference_model.data):
            # Increment the outer CV fold number after each iteration
            fold_number += 1
            print(f"{'=' * 50}\n{reference_model.regressor_name}, CV {fold_number}\n{'=' * 50}")
            # A new model object is initiated as a copy from the reference model
            model = deepcopy(reference_model)
            # Set config name to distinguish the different models in the output
            model.set_setup_name('CV{}'.format(fold_number))
            # Define the "outer loop" for a typical nested cross validation workflow
            model.define_outer_loop(train_index, test_index)
            # Gradually increase training fractions
            fractions_dict = model.get_train_data_fractions(
                np.arange(0.1, 1.1, 0.1))  # {0.1: [1,4,5], 0.2: [1,4,5,7,8], ...}
            for fraction in fractions_dict.keys():
                print(f'---------- Train fraction {fraction} -----------')
                this_model = deepcopy(model)
                # train on subset
                this_model.train_regressor_on_subset(subset_index=fractions_dict[fraction])
                # The optimized model is tested on the test set which was defined in the "outer loop"
                this_model.complete_evaluate_models(run_id=fold_number, train_fraction=fraction, visualize=False)

                self.collect_scores(this_model)

        v = Visualize(self, function_name)
        v.boxplot_scatterplot_CV_performance(by_categories=['train_fraction'])

    # -----------------------------------------------------#
    # Basic functions, duplicates without cross validation #
    # -----------------------------------------------------#
    def test_different_setups(self, regressor_list=None, feature_list=None,
                              split_by_compound=False, include_plant_fingerprints=False, number_of_splits=5,
                              load_existing=False, mode='singlevalue', config='default'):
        """
        Analogue to test_different_setups_cross_val() except that the split is random.
        Currently, we prefer to use cross validation instead of random splits.
        This split however might be important to test other tasks.
        For example see 'include_plant_fingerprints' parameter.

        :param regressor_list:
        :param feature_space_list:
        :param split_by_compound:
        :param include_plant_fingerprints:
        :param number_of_splits:
        :param load_existing:
        """
        print("\n############# Test different setups - NO nested cross validation #############")
        function_name = sys._getframe().f_code.co_name  # get the name of the function

        self.setting_the_stage(function_name, regressor_list, feature_list, mode='singlevalue', config='default',
                                                              load_existing=load_existing)
        # Iterate through feature spaces
        for feature_space in self.feature_space_list:
            self.descriptors.define_feature_space(feature_space)

            for regressor_dict in self.regressor_list:
                model = self.create_model_from_dict(regressor_dict)

                # This is where the method differs from other similar methods
                # -----------------------------------------------------------#
                if include_plant_fingerprints:
                    model.include_plant_fingerprints()
                model.simple_run(random_state_list=list(range(0, number_of_splits)),
                                 split_by_compound=split_by_compound,
                                 verbose=True)
                # -----------------------------------------------------------#

                self.collect_scores(model)

        self.visualize_test_scores(analysis_type=function_name, categories=['descriptors', 'regressor'])

    # -----------------------------------------------------#
    # Other functions and definitions                      #
    # -----------------------------------------------------#

    def setting_the_stage(self, function_name, regressor_name_list = None, feature_space_list=None,
                          mode='singlevalue', config='default', load_existing=False):
        """
        Preparing the modelling environment - includes setting output filenames, loading regressors from settings files,
        loading features
        :param function_name:
        :param regressor_name_list: List of keys as defined in regressor dictionaries, e.g., ['RF', 'SVR']
        :param feature_space_list: List of feature spaces to be explored, e.g., ['maccs', 'maccs+padel', 'all']
        :param mode: 'singlevalue' for defined settings or 'range' for grid search/optimization
        :param config: filename tag for loading regressor settings from config folder
        :param load_existing: load existing test scores for visualisation
        """
        print("-> setting the stage")
        self.clean_scores()
        self.build_scores_filenames(function_name)

        if load_existing:
            self.build_scores_filenames(function_name)
            self.test_scores = pd.read_csv(self.test_scores_tsv, sep='\t')
            # if regressor does not have scores, raise error
            if 'regressor_shortname' not in self.test_scores.columns.values:
                self.test_scores['regressor_shortname'] = Util.get_display_name(list(self.test_scores['regressor'].values))
            for reg in regressor_name_list:
                assert reg in self.test_scores['regressor_shortname'].values, f"No test scores found for {reg}"
            # if there are scores for a regressor not in regressor list, drop scores for visualization
            self.test_scores = self.test_scores[self.test_scores['regressor_shortname'].isin(regressor_name_list)]

        else:
            self.set_feature_space_list(feature_space_list)
            self.load_regressor_settings(mode=mode, config=config)
            self.set_regressor_list(regressor_name_list)
            if mode == 'range':
                self.expand_regressor_list_to_feature_selection()


    def prepare_reference_model(self):
        """
        Create a model from a regressor dictionary as defined in the config files
        :param regressor_dict: regressor dictionary
        :return: Model() object
        """
        # Define the pipeline
        reference_model = Model(self.pepper, self.descriptors, self.model_data)
        reference_model.prepare(verbose=True)
        return reference_model


    # ------- Regressors  related ----------- #

    def load_regressor_settings(self, mode ='singlevalue', config='default'):
        """
        Load regressor settings from yaml files from the config folder

        :param mode: 'singlevalue' or 'range' (for grid search)
        :param config: 'default', or user-defined regressor setting
        """
        filename = '../config/regressor_settings_{}_{}.yml'.format(mode, config)
        print("\tload regressor settings from {}".format(filename))
        with open(filename, 'r') as file:
            self.regressor_settings = yaml.safe_load(file)

        self.complete_regressor_name_list = list(self.regressor_settings.keys())

        # replace regressor string with function
        for reg in self.regressor_settings.keys():
            self.regressor_settings[reg]['regressor'] = self.get_regressor_by_name(self.regressor_settings[reg]['regressor'])




    def get_regressor_by_name(self, regressor_string):
        """
        Load regressor function from a regressor name
        :param regressor_string: name of regressor as defined in config file (function name with parentheses)
        :return: Regressor object
        """
        if regressor_string == 'RandomForestRegressor':
            return RandomForestRegressor(random_state=self.random_state)
        elif regressor_string == 'GradientBoostingRegressor':
            return GradientBoostingRegressor(random_state=self.random_state)
        elif regressor_string == 'AdaBoostRegressor':
            return AdaBoostRegressor(random_state=self.random_state)
        elif regressor_string == 'MLPRegressor':
            return MLPRegressor(random_state=self.random_state)
        elif regressor_string == 'SVR':
            return SVR()
        elif regressor_string == 'KNeighborsRegressor':
            return KNeighborsRegressor()
        elif regressor_string == 'GaussianProcessRegressor':
            return GaussianProcessRegressor(random_state=self.random_state)
        else:
            raise NotImplementedError('No regressor type defined for regressor_string = {}'.format(regressor_string))

    def set_regressor_list(self, regressor_list: list):
        """
        Sets a list of regressors for the modeling object. If called without providing a regressor list, then a default
        list (with all regressors currently implemented) will be assigned.
        :param regressor_list: e.g., ['RF', 'SVR']
        """
        self.regressor_name_list = regressor_list or self.complete_regressor_name_list
        for reg in self.regressor_name_list:
            assert self.regressor_settings.get(reg), 'Error: no regressor settings found for {}'.format(reg)
            self.regressor_list.append(self.regressor_settings.get(reg))

    def expand_regressor_list_to_feature_selection(self):
        """
        Load regressor collection and expands it where several feature selection methods are defined
        @param regressors: dictionary of regressors, as in self.regressor_name_list
        """
        self.regressor_list = []
        print('-> combinations to be tested:')
        for reg in self.regressor_name_list:
            for selection_method in self.regressor_settings[reg]["feature_selection_method"]:
                for reduction_method in self.regressor_settings[reg]["feature_reduction_method"]:
                    print('\t', reg, selection_method, reduction_method)
                    this_regressor = self.regressor_settings[reg].copy()

                    this_regressor["feature_selection_method"] = selection_method
                    if selection_method != "None": # ignore parameters
                        this_regressor["feature_selection_parameters"][selection_method] = \
                        self.regressor_settings[reg]["feature_selection_parameters"][selection_method]

                    this_regressor["feature_reduction_method"] = reduction_method
                    if reduction_method != "None": # ignore parameters
                        this_regressor["feature_reduction_parameters"][reduction_method] = \
                        self.regressor_settings[reg]["feature_reduction_parameters"][reduction_method]

                    self.regressor_list.append(this_regressor)

    def get_optimized_regressor_dict(self, config):
        """
        Load optimized regressor settings (singlevalue) from indicated config file and return a regressor dictionary
        :param config: user-defined regressor setting tag (e.g., 'soil_optimized')
        :return: dictionary with optimized regressor settings
        """
        self.load_regressor_settings(config=config)
        # ensure only one regressor is indicated for the optimized
        assert len(self.complete_regressor_name_list) == 1, (
            "Please provide only one regressor in the config file (now: {})".format(self.complete_regressor_name_list))
        # fetch regressor settings, feature space and feature selection and pass them to model
        regressor_dict = self.regressor_settings[self.complete_regressor_name_list[0]]
        return regressor_dict
        
    def get_regressor_name_list(self):
        return self.regressor_name_list

    # ------- Features related ----------- #

    def set_feature_space_list(self, feature_space_list):
        self.feature_space_list = feature_space_list or self.complete_feature_space_list

    # ------- Scores related ----------- #
    def clean_scores(self):
        self.test_scores = pd.DataFrame()
        self.train_scores = pd.DataFrame()

    def collect_scores(self, model: Model):
        """
        Collect scores from a Model object ad save it to the Modeling object

        @param model: model from where train_scores and test_scores are to be collected
        """
        # training
        self.train_scores = pd.concat([self.train_scores, model.train_scores_df], ignore_index=True)
        self.train_scores.to_csv(self.train_scores_tsv, sep='\t', index=True)
        # testing
        self.test_scores = pd.concat([self.test_scores, model.test_scores_df], ignore_index=True)
        self.test_scores.to_csv(self.test_scores_tsv, sep='\t', index=True)

    def build_scores_filenames(self, analysis_type: str):
        """
        Build filenames to save collected train and test scores

        @param analysis_type: name of analysis type, used to save file
        """
        # training and testing scores
        scores_dir = os.path.join(self.get_data_directory(), 'scores', analysis_type, self.curation_type)
        self.build_directory_structure(scores_dir)
        if not self.train_scores_tsv:
            self.train_scores_tsv = self.build_output_filename(os.path.join(scores_dir, 'train_scores'))
        if not self.test_scores_tsv:
            self.test_scores_tsv = self.build_output_filename(os.path.join(scores_dir, 'test_scores'))

    # ------- Visualizing ----------- #
    def visualize_test_scores(self, analysis_type: str, plot_type='boxplot_by_combination',
                              categories=['regressor', 'descriptors', 'feature_selection', 'feature_reduction']):
        """
        Visualise R2 and RMSE for combinations of regressors, descriptors, and/or feature selection methods
        :param analysis_type: name of the analysis type, used to save file
        :param plot_type: 'boxplot_by_combination' (default) or 'scatterplot_boxplot'
        (plots all combinations of 'categories')
        :param categories: list of categories for 'boxplot_by_combination'.
        Box plots are ordered and colored by categories[0].
        """
        v = Visualize(self, analysis_type)
        if plot_type == 'boxplot_by_combination':
            v.boxplot_CV_performance(categories=categories)
        elif plot_type == 'scatterplot_boxplot':
            v.boxplot_scatterplot_CV_performance()
        else:
            raise NotImplementedError('No plot type defined for plot_type = {}'.format(plot_type))


    def save_test_score_overview(self, overview_by= ['descriptors', 'regressor', 'combination']):
        """
        Save an overview of scores with averages per combination and per descriptor
        :param overview_by:
        """
        columns = overview_by + ['RMSE', 'R2']
        df = self.test_scores.loc[:,columns]
        for item in overview_by:
            self.save_scores_by(df, item)


    def save_scores_by(self, df, by):
        result_df = pd.DataFrame()
        result_df[by] = list(df.groupby(df[by])['R2'].mean().index)
        result_df['R2 mean'] = df.groupby(df[by])['R2'].mean().values
        result_df['R2 standard deviation'] = df.groupby(df[by])['R2'].std().values
        result_df['RMSE mean'] = df.groupby(df[by])['RMSE'].mean().values
        result_df['RMSE standard deviation'] = df.groupby(df[by])['RMSE'].std().values
        result_df = result_df.round(2)
        print_to = self.test_scores_tsv.replace('test_scores', f'test_scores_by_{by}')
        result_df.to_csv(print_to, sep='\t', index=False)
        print(f'-> Save scores summarized by {by} to {print_to}')