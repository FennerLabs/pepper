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
from sklearn.model_selection import KFold, train_test_split
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
import seaborn as sns
import matplotlib.pyplot as plt


# Importing all regressors
# from sklearn.linear_model import LinearRegression
# from sklearn.linear_model import Ridge
# from sklearn.linear_model import SGDRegressor
# from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel
from sklearn.neural_network import MLPRegressor

# from sklearn.tree import DecisionTreeRegressor


class Modeling(Pepper):
    def __init__(self, pep: Pepper(), data: DataStructure, descriptors: Descriptors):
        super().__init__()
        self.set_data_directory(os.path.join(pep.data_directory, 'modeling'))
        self.tag = pep.get_tag()
        self.setup_name = pep.get_setup_name()
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

        # Scoring storage data objects
        self.test_scores = pd.DataFrame()
        self.train_scores = pd.DataFrame()
        self.holdout_scores = pd.DataFrame()

        # filenames
        self.train_scores_tsv = None
        self.test_scores_tsv = None
        self.holdout_scores_tsv = None

        # nested folders with scores by model setup
        self.train_scores_tsv_dict = {}
        self.test_scores_tsv_dict = {}
        self.holdout_scores_tsv_dict = {}

        self.predicted_target_variable_train_tsv_dict = {}
        self.predicted_target_variable_test_tsv_dict = {}
        self.predicted_target_variable_holdout_tsv_dict = {}

        self.predicted_target_variable_tsv = None

        # cross validation predictions
        self.predicted_target_variable_train = pd.DataFrame()
        self.predicted_target_variable_test = pd.DataFrame()
        self.predicted_target_variable_train_tsv = None
        self.predicted_target_variable_test_tsv = None

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

    def run_test_model(self, regressor_name = 'RF', feature_space = 'maccs'):
        """
        # Run a test model to verify correct data and descriptor config. Used for teaching purposes mainly.
        """
        print("############# Run test model #############")
        # get the name of the function
        function_name = sys._getframe().f_code.co_name
        # define which descriptors to use
        self.descriptors.define_feature_space(feature_space)
        # setting the stage
        self.setting_the_stage(function_name=function_name, regressor_name_list=[regressor_name],
                               feature_space_list=[feature_space], mode='singlevalue', config='default',
                               load_existing=False)
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

    def cross_validation_screening(self, number_of_splits=5,
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

        if not load_existing:
            # Iterating through feature spaces to be explored
            assert len(feature_space_list) > 0, 'Feature space list cannot be empty'
            for feature_space in self.feature_space_list:
                self.descriptors.define_feature_space(feature_space)
                reference_model = self.prepare_reference_model()
                # Create a KFold object with n_splits
                kf = KFold(n_splits=number_of_splits, shuffle=True, random_state=self.random_state)
                # Iterate through regressor list
                for regressor_dict in self.regressor_list:
                    reference_model.load_settings_from_dict(regressor_dict)
                    reference_model.settings_string += f'_{feature_space}'
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
            self.write_scores_filesystem()
            # print test score overview
            self.save_test_score_overview(overview_by=['descriptors','regressor'])
        self.visualize_test_scores(analysis_type=function_name, categories=['descriptors','regressor'])  #



    def nested_cross_validation_screening(self, regressor_name_list = ['RF', 'GB','AB', 'SVR', 'GPR', 'KNN'],
                                   feature_space_list = ['all'],
                                   load_existing=False, config='default'):
        """
        This method intends to easily evaluate many combinations of regressors and descriptors
        to understand which (features,regressor) pairs perform better for a given dataset.
        It is done in a typical nested cross validation fashion.
        Both "inner" and "outer" loops are 5-fold-CV

        :param regressor_name_list: List of abbreviations of regressors for which nested CV will be performed
        :param feature_space_list: Features to be considered. By default, all loaded features are considered
        :param load_existing: Load existing score files or visualisation
        :param config: regressor configs as defined in ./config/regressor_settings__['range'/'singlevalue']_[config].yml config file
        """
        function_name = sys._getframe().f_code.co_name  # get the name of the function
        print("\n############# Nested cross-validation screening #############")
        self.setting_the_stage(function_name, regressor_name_list=regressor_name_list, feature_space_list=feature_space_list,
                               mode='range', config=config, load_existing=load_existing)

        if not load_existing:
            # Iterating through feature spaces to be explored
            assert len(feature_space_list) > 0, 'Feature space list cannot be empty'
            for feature_space in feature_space_list:
                self.descriptors.define_feature_space(feature_space)

                # Create reference Model object, prepare data and feature matrices, preprocess features
                reference_model = self.prepare_reference_model()

                # split of a 10% holdout set to be tested but only split the indices
                # self.train_data, self.holdout_data = train_test_split(reference_model.data, test_size=0.1, random_state=self.random_state)
                # Create a KFold object with n_splits
                kf = KFold(n_splits=5, shuffle=True, random_state=self.random_state)

                # perform nested CV with feature selection, hyperparameter tuning and evaluation on test set for each
                # combination of regressor and feature selection method
                for regressor_dict in self.regressor_list:
                    # Get the settings of the regressor
                    reference_model.load_settings_from_dict(regressor_dict)
                    # indicate feature_space in settings string
                    reference_model.settings_string += f'_{feature_space}'
                    # clear best models list
                    self.best_models = []
                    # Nested cross validation "outer loop"
                    fold_number = 0
                    feature_reduction = regressor_dict['feature_reduction_method']

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
                        model.complete_train_regressors(run_id=fold_number)

                        # The optimized model is tested on the test set which was defined in the "outer loop"
                        model.complete_evaluate_models(run_id=fold_number, get_nearest_neighbors=True)

                        # Best settings are saved and performance scores collected
                        self.best_models.append((model.regressor_name, model.best_regressor_model, model.best_reduction_params))
                        self.collect_scores(model)

                    # # train the best model on the entire training set and evaluate on holdout set with best parameters
                    # model = self.train_and_evaluate_on_holdout_set(reference_model, regressor_dict)

                    # self.collect_scores_holdout(model, feature_space=feature_space)

                    # save the predicted target variable for all cross validation runs
                    self.save_predictions_cross_validation(model, feature_reduction, feature_space=feature_space)

                    self.visualize_predictions(analysis_type=function_name, prediction_type='cv', model=model, feature_space=feature_space, feature_reduction=feature_reduction)
                    # self.visualize_predictions(analysis_type=function_name, prediction_type='holdout', model=model, feature_space=feature_space, feature_reduction=feature_reduction)

                    # Clear dataframes for next run

                    self.predicted_target_variable_train = pd.DataFrame()
                    self.predicted_target_variable_test = pd.DataFrame()
                    self.predicted_target_variable_holdout = pd.DataFrame()
            self.write_scores_filesystem()
            # print test score overview
            self.save_test_score_overview(['regressor', 'descriptors'])

        # Visualisation of performance scores of the optimized models on the respective test sets
        self.visualize_test_scores(analysis_type=function_name,
                                   categories=['regressor', 'descriptors', 'feature_selection', 'feature_reduction'])


    def collect_scores_holdout(self, model, feature_space: str = 'all'):
        """
        Collect holdout scores after evaluation on holdout set
        """
        self.holdout_scores =  pd.concat([self.holdout_scores, model.holdout_scores_df], ignore_index=True)
        self.holdout_scores.to_csv(self.holdout_scores_tsv_dict[model.regressor_name_short][feature_space], sep='\t', index=True)
        self.predicted_target_variable_holdout = model.predicted_target_variable_holdout
    
    def save_predictions_cross_validation(self, model, feature_reduction, feature_space: str = 'all'):
        """
        Save the predicted target variable for all cross validation runs
        """
        
        
        if not self.predicted_target_variable_train.empty:
            self.predicted_target_variable_train.to_csv(self.predicted_target_variable_train_tsv_dict[model.regressor_name_short][feature_space][feature_reduction], sep='\t', index=False)

        if not self.predicted_target_variable_test.empty:
            self.predicted_target_variable_test.to_csv(self.predicted_target_variable_test_tsv_dict[model.regressor_name_short][feature_space][feature_reduction], sep='\t', index=False)

        if hasattr(self, "predicted_target_variable_holdout") and not self.predicted_target_variable_holdout.empty:
            self.predicted_target_variable_holdout.to_csv(self.predicted_target_variable_holdout_tsv_dict[model.regressor_name_short][feature_space][feature_reduction], sep='\t', index=False)

    # -----------------------------------------------------#
    # Functions for optimized regressors  #
    # -----------------------------------------------------#
    def build_gpytorch_model(self, regressor_name: str, feature_space: str, config):
        function_name = sys._getframe().f_code.co_name  # get the name of the function
        self.setting_the_stage(function_name, regressor_name_list = [regressor_name],
                               feature_space_list=[feature_space], mode='singlevalue', config=config)
        regressor_dict = self.get_optimized_regressor_dict(config, regressor_name)
        
        model = Model(self.pepper, self.descriptors, self.model_data)
        model.load_settings_from_dict(regressor_dict)
        model.descriptors.define_feature_space(feature_space)
        model.prepare()
        model.use_all_data_for_training()

        model.select_features()
        model.reduce_feature_space()
        model.set_regressor_parameters(model.regressor, model.regressor_params)
        return model

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
        regressor_dict = self.get_optimized_regressor_dict(config, regressor_name)
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

        model.save_model()

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
                               feature_space_list=[feature_space], mode='singlevalue', config=config,
                               load_existing=load_existing)

        if not load_existing:
            regressor_dict = self.get_optimized_regressor_dict(config, regressor_name)
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
        print("-> setting the stage:")
        self.clean_scores()

        if load_existing:
            print('\tloading existing scores for visualization')
            self.build_scores_filenames(function_name)
            self.load_scores()
            # if regressor does not have scores, raise error
            if 'regressor_shortname' not in self.test_scores.columns.values:
                self.test_scores['regressor_shortname'] = Util.get_display_name(list(self.test_scores['regressor'].values))
            for reg in regressor_name_list:
                assert reg in self.test_scores['regressor_shortname'].values, f"No test scores found for {reg}"
            # if there are scores for a regressor not in regressor list, drop scores for visualization
            self.test_scores = self.test_scores[self.test_scores['regressor_shortname'].isin(regressor_name_list)]

        else:
            print('\tsetting model configurations')
            self.set_feature_space_list(feature_space_list)
            self.load_regressor_settings(mode=mode, config=config)
            self.set_regressor_list(regressor_name_list)
            self.build_scores_filenames(function_name)
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

        # Ensure the settings are defined correctly:
        mandatory_settings = ['name', 'regressor', 'regressor_params']
        optional_settings = ['feature_reduction_method', 'feature_selection_method', 'feature_reduction_parameters', 'feature_selection_parameters']
        for regressor, setting in self.regressor_settings.items():
            for mandatory in mandatory_settings:
                assert mandatory in setting.keys(), (f'Mandatory setting {mandatory} is missing '
                                                     f'from {regressor} regressor in config file')
            for option in optional_settings:
                if option not in setting.keys():
                    if mode == 'singlevalue':
                        self.regressor_settings[regressor][option] = 'None'
                    else:
                        self.regressor_settings[regressor][option] = ['None']

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
            return GaussianProcessRegressor(kernel=ConstantKernel(2.0, (1e-3, 1e3)) * Matern(length_scale=2.5, length_scale_bounds=(1e-3, 1e3)), random_state=self.random_state)
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
        self.regressor_list = [] # re-initialize regressor list
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

    def get_optimized_regressor_dict(self, config, regressor_name):
        """
        Load optimized regressor settings (singlevalue) from indicated config file and return a regressor dictionary
        :param config: user-defined regressor setting tag (e.g., 'soil_optimized')
        :param config: regressor name (to select parameters from config file)
        :return: dictionary with optimized regressor settings
        """
        self.load_regressor_settings(config=config)
        # ensure only one regressor is indicated for the optimized
        assert regressor_name in self.complete_regressor_name_list, (
            "Please provide settings for the {regressor_name} regressor in the config file")
        # fetch regressor settings, feature space and feature selection and pass them to model
        regressor_dict = self.regressor_settings[regressor_name]
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
        ### global scores files ###
        # training
        self.train_scores = pd.concat([self.train_scores, model.train_scores_df], ignore_index=True)
        self.train_scores.to_csv(self.train_scores_tsv, sep='\t', index=True)
        self.predicted_target_variable_train = pd.concat([self.predicted_target_variable_train, model.predicted_target_variable_train], ignore_index=True)
        self.predicted_target_variable_train.to_csv(self.predicted_target_variable_train_tsv, sep='\t', index=True)

        # testing
        self.test_scores = pd.concat([self.test_scores, model.test_scores_df], ignore_index=True)
        self.test_scores.to_csv(self.test_scores_tsv, sep='\t', index=True)
        self.predicted_target_variable_test = pd.concat(
            [self.predicted_target_variable_test, model.predicted_target_variable_test], ignore_index=True)
        self.predicted_target_variable_test.to_csv(self.predicted_target_variable_test_tsv, sep='\t', index=True)

    def write_scores_filesystem(self):
        for regressor in self.regressor_list:
            regressor_name = regressor["name"]
            for feature_space in self.feature_space_list:
                df_test_regressors =  self.test_scores[self.test_scores['regressor'] == regressor_name]
                df_train_regressors = self.train_scores[self.train_scores['regressor'] == regressor_name]
                df_test = df_test_regressors[df_test_regressors['descriptors'] == feature_space]
                df_train = df_train_regressors[df_train_regressors['descriptors'] == feature_space]
                df_test.to_csv(self.test_scores_tsv_dict[Util.get_display_name(regressor_name)][feature_space], sep='\t', index=True)
                df_train.to_csv(self.train_scores_tsv_dict[Util.get_display_name(regressor_name)][feature_space], sep='\t', index=True)




    def build_scores_filenames(self, analysis_type: str):
        """
        Build filenames to save collected train and test scores

        @param analysis_type: name of analysis type, used to save file
        """
        # training and testing scores
        scores_dir = os.path.join(self.get_data_directory(), 'scores', analysis_type, self.curation_type)
        self.build_directory_structure(scores_dir)

        # for overview of scores
        self.test_scores_tsv = os.path.join(scores_dir, f'test_scores.tsv')
        self.train_scores_tsv = os.path.join(scores_dir, f'train_scores.tsv')

        # create subfolders for each regressor
        score_types = {
            'train_scores': self.train_scores_tsv_dict,
            'test_scores': self.test_scores_tsv_dict,
            'holdout_scores': self.holdout_scores_tsv_dict,
            'predicted_target_variable_train': self.predicted_target_variable_train_tsv_dict,
            'predicted_target_variable_test': self.predicted_target_variable_test_tsv_dict,
            'predicted_target_variable_holdout': self.predicted_target_variable_holdout_tsv_dict,
        }

        # Loop through regressors and feature spaces and feature reduction methods to build filenames
        for regressor in self.regressor_name_list:
            regressor_dir = os.path.join(scores_dir, regressor)
            self.build_directory_structure(regressor_dir)

            for features_space in self.feature_space_list:
                for score_name, score_dict in score_types.items():
                    # Ensure nested dicts are initialized
                    score_dict.setdefault(regressor, {})

                    # only add feature reduction method for predicted scores
                    if 'predicted' in score_name:
                        # Ensure nested dicts are initialized
                        score_dict[regressor].setdefault(features_space, {})
                        feature_reduction_methods = self.regressor_settings[regressor]['feature_reduction_method']

                        if type(feature_reduction_methods) == str:
                            feature_reduction_methods = [feature_reduction_methods]

                        for fr in feature_reduction_methods:
                            path = os.path.join(regressor_dir, f"{score_name}_{features_space}_{fr}_{self.get_setup_name()}")
                            score_dict[regressor][features_space][fr] = self.build_output_filename(path)
                    else:
                        path = os.path.join(regressor_dir, f"{score_name}_{features_space}_{self.get_setup_name()}")
                        score_dict[regressor][features_space] = self.build_output_filename(path)

    def load_scores(self):
        # Scores are loaded from the overall train and test scores files (not the nested ones)
        if self.test_scores_tsv:
            self.test_scores = pd.read_csv(self.test_scores_tsv, sep='\t')
        if self.train_scores_tsv:
            self.train_scores = pd.read_csv(self.train_scores_tsv, sep='\t')


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



    def visualize_predictions(self, analysis_type: str, prediction_type='cv', model=None, feature_space: str = 'all', feature_reduction: str = 'None'):
        v = Visualize(self, analysis_type)
        if prediction_type == 'cv':
            v.parity_predicted_vs_test(prediction_type='cv', model=model, feature_space=feature_space, feature_reduction=feature_reduction)
            if model.has_y_pred_score:
                v.calibration_plots(prediction_type='cv', model=model, feature_space=feature_space, feature_reduction=feature_reduction)
                v.uncertainty_distribution_plots(prediction_type='cv', model=model, feature_space=feature_space, feature_reduction=feature_reduction)
        elif prediction_type == 'holdout':
            v.parity_predicted_vs_test(prediction_type='holdout', model=model, feature_space=feature_space, feature_reduction=feature_reduction)
            if model.has_y_pred_score:
                v.calibration_plots(prediction_type='holdout', model=model, feature_space=feature_space, feature_reduction=feature_reduction)
                v.uncertainty_distribution_plots(prediction_type='holdout', model=model, feature_space=feature_space, feature_reduction=feature_reduction )

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

                    # a function to get the best parameters over all 5 folds
    def get_best_parameters(self, reference_model):
        # Collect best parameters from each fold
        best_params_list = [model[1].get_params() for model in self.best_models if model[0] == reference_model.regressor_name]

        best_reduction_params_list = [{"reduction_n": model[2]} for model in self.best_models if model[0] == reference_model.regressor_name and len(model) > 1]
        # filter out array like values

        # append best reduction params to best params
        for p in range(len(best_params_list)):
            best_params_list[p].update(best_reduction_params_list[p])

        clean_params_list = []
        for p in best_params_list:
            clean = {k: v for k, v in p.items() if not isinstance(v, (list, np.ndarray))}
            
            clean_params_list.append(clean)

        # Calculate the frequency of each parameter set
        params_df = pd.DataFrame(clean_params_list)
       

        most_common_params = params_df.mode().iloc[0].to_dict()
        for col in params_df.columns:
            try:
                most_common_params[col] = params_df[col].mode(dropna=True).iloc[0]
            except Exception:
                most_common_params[col] = params_df[col].dropna().iloc[0] if not params_df[col].dropna().empty else None

        # convert max_depth to int if exists
        if 'max_depth' in most_common_params:
            most_common_params['max_depth'] = int(most_common_params['max_depth'])
        return most_common_params
    
    def train_and_evaluate_on_holdout_set(self, reference_model, regressor_dict):
        # A new model object is initiated as a copy from the reference model
        model = deepcopy(reference_model)
        model.set_setup_name('holdout')
        model.train_indices = self.train_data.index
        model.test_indices = self.holdout_data.index
        model.X_train = model.features.iloc[model.train_indices]
        model.X_test = model.features.iloc[model.test_indices]
        model.smiles_train = model.data[model.smiles_name].iloc[model.train_indices]
        model.smiles_test = model.data[model.smiles_name].iloc[model.test_indices]
        model.y_train = model.target_variable.iloc[model.train_indices][self.target_variable_name]
        model.y_test = model.target_variable.iloc[model.test_indices][self.target_variable_name]
        if self.has_target_variable_std():
            model.y_train_std = model.target_variable_std.iloc[model.train_indices]
            model.y_test_std = model.target_variable_std.iloc[model.test_indices]

        best_params = self.get_best_parameters(reference_model)
        # find the regressor parameters set in the regressor dict
        model.best_regressor_params = {}
        for key, value in regressor_dict['regressor_params'].items():
            if key in best_params.keys():
                model.best_regressor_params[key] = best_params[key]
        model.set_regressor_parameters(model.regressor, model.best_regressor_params)
        model.no_feature_selection()
        #TODO feature reduction
        if model.feature_reduction_method not in [None, 'None']:
            model.reduce_feature_dimensionality(n_components=best_params.get('reduction_n'), method_name='PCA')
            model.best_reduction_params = best_params.get('reduction_n')

        else:
            model.no_feature_reduction()
        model.train_model()
        model.predict(model.X_test_reduced,dataset_type='test')
        # model.calculate_holdout_scores()
        model.save_scores(stage='hold_out')
        model.save_predicted_values(stage='hold_out', random_state=self.random_state, get_nearest_neighbors=True)
        return model
