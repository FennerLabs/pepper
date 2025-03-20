import os

import pandas as pd
import numpy as np
import random
import sys
from copy import deepcopy

from mlxtend.preprocessing import standardize
from mordred import Descriptor
from pandas.core.common import random_state
from sklearn.model_selection import train_test_split, cross_validate
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
# from sklearn.feature_selection import SequentialFeatureSelector # this library does not keep CV scores
from mlxtend.feature_selection import SequentialFeatureSelector
from sklearn.ensemble import RandomForestRegressor
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD

from genetic_selection import GeneticSelectionCV

from sklearn.preprocessing import MinMaxScaler, FunctionTransformer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import VarianceThreshold
# from sklearn.feature_selection import SelectKBest
# from sklearn.feature_selection import mutual_info_regression

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from collections import defaultdict

from pepper_lab.pepper import Pepper
from pepper_lab.visualize import Visualize


class Model(Pepper):
    def __init__(self, pep: Pepper, descriptors, model_data, pipe=None, regressor=None):
        """ features: modeling.features
        target_variable: modeling.target_variable"""
        super().__init__()

        self.tag = pep.get_tag()
        self.set_data_directory(os.path.join(pep.data_directory, 'modeling', 'models'))
        self.section = 'model'
        self.data_type = pep.get_data_type()
        self.setup_name = pep.get_setup_name()
        self.target_variable_name = pep.get_target_variable_name()
        self.target_variable_std_name = pep.get_target_variable_std_name()
        self.compound_name = pep.get_compound_name()
        self.id_name = pep.get_id_name()
        self.curation_type = pep.get_curation_type()
        self.pepper = pep
        self.random_state = pep.get_random_state()

        # Regressor and pipes
        self.regressor = regressor
        self.pipe = pipe
        self.regressor_params = {}
        self.best_regressor_params = {}
        self.best_regressor_model = None
        self.test_size = 0.2
        # feature selection method and parameters (e.g., sequential feature selection, genetic algorithm)
        self.feature_selection_method = None
        self.feature_selection_params = {}
        # feature dimensionality reduction methods and parameters (e.g., PCA: Principal Component Analysis)
        self.feature_reduction_method = None
        self.feature_reduction_params = {}

        # data obtained from input
        self.descriptors = descriptors
        self.data = model_data
        self.target_variable = model_data[self.target_variable_name]
        self.smiles_name = descriptors.smiles_name

        # modeling variables, set by self.prepare()
        self.features = pd.DataFrame()  # This could be modeling.features or modeling.reduced_features

        # data variables
        self.X = pd.DataFrame()  # 'defined after applying get_X_y' #todo: remove, not used?
        self.y = pd.DataFrame()  # 'defined after applying get_X_y' #todo: remove, not used?
        self.X_train = pd.DataFrame()  # 'X_train not defined yet; apply data_split'
        self.X_test = pd.DataFrame()  # 'X_test not defined yet; apply data_split'
        self.y_train = ()  # 'y_train not defined yet; apply data_split'
        self.y_test = ()  # 'y_test not defined yet; apply data_split'
        self.y_train_std = ()  # standard deviation of y obtained from BI, uncertainty of y_train
        self.y_test_std = () # standard deviation of y obtained from BI, uncertainty of y_test
        self.test_indices = [] # indices of self.data used for testing
        self.train_indices = []  # indices of self.data used for training

        # Feature reduction and selection
        self.feature_preprocessor = None # use to preprocess X - self.feature_preprocessor.transform(X)
        self.feature_transformer = None # use to obtain X_selected from X - self.feature_transformer.transform(X)
        self.feature_names = []
        # X after feature selection
        self.X_train_selected = pd.DataFrame()
        self.X_test_selected = pd.DataFrame()
        # X after feature dimensionality reduction
        self.X_train_reduced = pd.DataFrame()
        self.X_test_reduced = pd.DataFrame()


        # To export results
        # original data as loaded from source
        self.y_pred_train = pd.DataFrame() # 'model not trained yet or predict method not applied to X_train'
        self.y_pred = pd.DataFrame() #'predict method not applied yet to X_test'
        self.y_pred_score = [] # predicted score for each value of y, e.g., stdev in case of GPR
        self.has_y_pred_score = False # True if a prediction score is expected from the regressor
        # fitting values
        self.fitting_data_values = pd.DataFrame()
        self.fitting_data_tsv = self.build_output_filename('training')
        # predicted values
        self.predicted_target_variable = pd.DataFrame()
        self.predicted_target_variable_tsv = self.build_output_filename('predictions')
        self.training_dict = {}

        # scores
        self.train_scores = pd.DataFrame()
        self.test_scores = pd.DataFrame()

        # To save as dataframes and to export as tsv
        self.train_scores_df = pd.DataFrame()
        self.train_scores_df_tsv = self.build_output_filename('train_scores')
        self.test_scores_df = pd.DataFrame()
        self.test_scores_df_tsv = self.build_output_filename('test_scores')

        # To save optimized settings and selected features
        self.optimized_settings_txt = ''
        self.selected_features_scores = pd.DataFrame(columns=['Feature', 'Score', 'CV_run', 'Average_score'])
        self.selected_features_scores_tsv = ''
        self.selected_features_names = []
        self.settings_string = '' # regressor name and feature selection method used to build model
        self.feature_names_used_for_training = [] # list of features needed to apply model to external data

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

    ################################## prepare model setup ###################################
    def build_output_filename(self, filename_string, suffix='.tsv'):
        """
        Create filename based on current directory, data typ, and tag
        :param filename_string: file name as a string
        :param suffix: file type, default is tab separated values
        :return:
        """
        file_string = (filename_string + '_{}_{}_{}' + suffix).format(self.data_type, self.tag, self.setup_name)
        complete_filename = os.path.join(self.data_directory, file_string)
        return complete_filename

    def prepare(self, verbose=True):
        """
        This function sets the regressor, the data and the features of the Model object.
        It also preprocesses the features (dropping NAs, scaling, removing highly correlated features)
        @param verbose: print preparation steps
        """
        if verbose:
            print('-> prepare model')

        # check regressor settings and put defaults if necessary
        if self.regressor is None: # set default
            print('Warning: using default regressor RF')
            self.regressor = RandomForestRegressor(random_state=self.random_state)
            self.regressor_name = 'RandomForestRegressor'

        self.features = self.descriptors.features
        # remove feature columns that have NA for one or more compounds
        self.remove_na_features(verbose=verbose)

        # get list of SMILESs that are common between data and descriptors
        smiles_list = list(self.features[self.smiles_name])

        # ensure that self.features and self.data contain the same set of compounds
        self.data = self.data[self.data[self.smiles_name].isin(smiles_list)]
        self.features = pd.merge(self.data[[self.smiles_name]], self.features, on=self.smiles_name, how='left')

        # data and features are re-indexed to avoid invalid indexes when splitting using KFold
        self.data.reset_index(drop=True, inplace=True)
        self.features.reset_index(drop=True, inplace=True)

        # set target variable from self.data
        self.target_variable = self.data[[self.target_variable_name]]
        # also set target_variable_std if a standard deviation is given
        if self.has_target_variable_std():
            self.target_variable_std = self.data[self.target_variable_std_name]
        else:
            self.target_variable_std = len(self.target_variable) * [np.NaN]

        # Drop SMILES column here
        self.features.drop(columns=self.smiles_name, inplace=True)

        # preprocess and preselect features
        self.preprocess_features()   # Here features are scaled: self.features -> self.features
        self.preselect_features()   # Here highly correlated and low variance are removed: self.features -> features


    def include_plant_fingerprints(self):
        self.features = pd.merge(self.descriptors.plant_fp, self.features, how='left',
                                 left_index=True, right_index=True)
        self.features.drop(columns=['plant'], inplace=True)
        print("Shape including the fingerprints: {}".format(self.features.shape))

    def define_regressor_by_dict(self, regressor_dict: dict):
        """
        Define regressor by dictionary
        :param regressor_dict: Dictionary in the format {'regressor': regressor object, 'name': regressor name,
        'regressor_params': {'param1': [1,2,3], 'param2':[True]}, 'feature_selection_method': 'importance',
        'feature_selection_params': {'method': 'top_k', 'top_k': 5}, 'descriptors': 'all'}
        """
        # define regressor settings
        self.regressor = regressor_dict['regressor']
        if regressor_dict.get('name'):
            self.regressor_name = regressor_dict['name']
        else:
            self.regressor_name = str(regressor_dict['name']).split('(')[0]
        if regressor_dict.get('regressor_params'):
            self.regressor_params = regressor_dict['regressor_params']
        # define feature space
        if regressor_dict.get('descriptors'):
            self.descriptors.define_feature_space(regressor_dict['descriptors'])
        # Define feature selection method and parameters
        if regressor_dict.get('feature_selection_method'):
            self.feature_selection_method = regressor_dict['feature_selection_method']
            if regressor_dict.get('feature_selection_parameters') and self.feature_selection_method != 'None':
                self.feature_selection_params = regressor_dict['feature_selection_parameters'][self.feature_selection_method]
        # Define feature reduction method and parameters
        if regressor_dict.get('feature_reduction_method'):
            self.feature_reduction_method = regressor_dict['feature_reduction_method']
            if regressor_dict.get('feature_reduction_parameters') and self.feature_reduction_method != 'None':
                self.feature_reduction_params = regressor_dict['feature_reduction_parameters'][self.feature_reduction_method]
        # define if output score is expected for this regressor
        if self.regressor_name in ['Gaussian Process Regressor', 'KNN Regressor']:
            self.has_y_pred_score =True

    def load_settings_from_dict(self, regressor_dict):
        self.define_regressor_by_dict(regressor_dict)
        self.set_feature_preselection_pipe()

    ############### Feature preprocessing, pre-selection, selection, transformation #########################

    def remove_na_features(self, verbose=True):
        """
        Removes the features that could not be calculated for ALL compounds.
        This is the simplest way to remove NA values but different methods
        to fill these values (instead of simply dropping) should be considered.
        """
        # Save smiles for putting them back later
        SMILES_series = self.features[[self.smiles_name]]

        # Consider possible "strange" values and drop
        if verbose:
            print("-> Removing features that were not calculated for all compounds")
            print("\tn features before: {}".format(len(self.features.columns)-1)) # don't count the smiles column
        self.features = self.features.apply(pd.to_numeric, errors='coerce')
        self.features.replace(r'', np.nan, regex=True)  # empty strings to nan #Question @Jose: why a regex here?
        self.features.replace([np.inf, -np.inf], np.nan, inplace=True)  # infinite to nan
        self.features.dropna(axis='index', how='all', inplace=True)
        self.features.dropna(axis='columns', how='any', inplace=True)

        # Place back SMILES
        self.features = self.features.merge(SMILES_series, how='right', right_index=True, left_index=True)

        if verbose:
            print("\tn features after: {}".format(len(self.features.columns)-1))  # don't count the smiles column

    def preprocess_features(self):
        """
        Min-max scaling of features
        """
        print("-> preprocess features (min-max scaling)")
        self.feature_preprocessor = MinMaxScaler().set_output(transform='pandas')
        self.feature_preprocessor.fit(self.features)
        self.features = self.feature_preprocessor.transform(self.features)

    def preselect_features(self):
        """
        Remove features with low variance and highly correlated features
        """
        print("-> pre-select features")
        print("\tn features before pre-selection: {}".format(len(self.features.columns)))
        self.set_feature_preselection_pipe()
        self.feature_preselector.fit(self.features)
        self.features = self.feature_preselector.transform(self.features)
        self.features_names = self.features.columns
        print("\tn features after pre-selection: {}".format(len(self.features.columns)))

    def set_feature_preselection_pipe(self):
        self.feature_preselector = Pipeline([
            ('variance_threshold', VarianceThreshold(threshold=0.02).set_output(transform='pandas')),
            # Feature selection based on variance threshold
            ('custom_feature_selection', FunctionTransformer(func=self.remove_highly_correlated_features,
                                                             validate=False,
                                                             kw_args={'corr_method': 'spearman',
                                                                      'cluster_threshold': 0.01}))
        ])

    @staticmethod
    def remove_highly_correlated_features(X_train, corr_method: str = 'spearman',
                                          cluster_threshold: float = 0.01, ignore=False):
        if ignore:
            return X_train
            # pass
        else:
            print('\t-> applying cluster threshold of {}'.format(cluster_threshold))
            print('\t\tn features before: {}'.format(len(X_train.columns)))
            corr = X_train.corr(corr_method).to_numpy()

            # Ensure the correlation matrix is symmetric
            corr = (corr + corr.T) / 2
            np.fill_diagonal(corr, 1)
            corr = np.nan_to_num(corr)

            # code from https://scikit-learn.org/stable/auto_examples/inspection/
            # plot_permutation_importance_multicollinear.html
            # We convert the correlation matrix to a distance matrix before performing
            # hierarchical clustering using Ward's linkage.
            distance_matrix = 1 - np.abs(corr)
            dist_linkage = hierarchy.ward(squareform(distance_matrix))

            cluster_ids = hierarchy.fcluster(dist_linkage, cluster_threshold, criterion="distance")
            cluster_id_to_feature_ids = defaultdict(list)
            for idx, cluster_id in enumerate(cluster_ids):
                cluster_id_to_feature_ids[cluster_id].append(idx)
            my_selected_features = [v[0] for v in cluster_id_to_feature_ids.values()]
            my_selected_features_names = X_train.columns[my_selected_features]
            X_train_sel = X_train[my_selected_features_names]
            print('\t\tn features after: {}'.format(len(my_selected_features_names)))
            return X_train_sel

    def select_features(self):
        kwargs = self.feature_selection_params
        self.settings_string = f"{self.regressor_name}_{self.feature_selection_method}"
        if self.feature_selection_method == 'importance':
            self.select_features_by_importance(**kwargs)
        elif self.feature_selection_method == 'sequential':
            self.select_features_by_sequential_selector(**kwargs)
        elif self.feature_selection_method == 'genetic':
            self.select_features_by_genetic_algorithm(**kwargs)
        elif self.feature_selection_method in [None, 'None']:
            self.no_feature_selection()
        else:
            raise ValueError(
                "Possible values for feature_selection parameter: 'importance', 'genetic', 'sequential', None")

        print(f"\tselected features: {len(self.selected_features_names)}")
        print(self.selected_features_names)
        # features that need to be calculated for endpoint prediction are saved here
        self.feature_names_used_for_training = self.selected_features_names

    def no_feature_selection(self):
        print(f'-> No feature selection applied')
        self.X_train_selected = self.X_train
        self.X_test_selected = self.X_test
        self.selected_features_names = self.X_train.columns

    def select_features_by_importance(self, method ='importance', top_k = 50, min_importance = 0.001, **kwargs):
        """
        Select features by importance from ensemble methods
        @param method: 'top_k' (default) or 'importance'
        @param top_k: Number of features to be selected. Set to 40 for 'top_k' default.
        @param min_importance: Minimal importance score to select feature. Set to any value between 0 and 1.
        """

        print(f'-> Select features by importance. Method: { method }')

        # Compute the feature importances after preprocessing
        cv_results = cross_validate(self.regressor, self.X_train, self.y_train, cv=5, scoring='r2', return_estimator=True)

        # save mean feature scores for feature selection
        feature_importances = np.mean(
            [est.feature_importances_ for est in cv_results['estimator']], axis=0)

        # Create a dictionary mapping feature importances to feature names
        feature_importance_map = {feature_name: importance for feature_name, importance in
                                  zip(self.features_names, feature_importances)}

        # Select features
        if method == 'top_k':
            # first, check if enough features available for feature selection
            if not self.has_enough_features_for_selection(top_k):
                return
            print(f'-> Selecting top {top_k} features')
            selected_features = np.argsort(feature_importances)[::-1][:top_k]
            self.selected_features_names = self.features_names[selected_features].values
            self.settings_string += f"_k{top_k}"
        elif method == 'importance':
            print(f'-> Selecting features with importance > {min_importance}')
            for k in feature_importance_map.keys():
                if feature_importance_map[k] >= min_importance:
                    self.selected_features_names.append(k)
            self.settings_string += f"_min{min_importance}"
        else:
            raise ValueError("The 'method' parameter should be set to 'top_k' or 'importance'")

        # save cv importance scores for visualization
        for index, cv_run in enumerate(cv_results['estimator']):
            df_run = pd.DataFrame()
            df_run['Feature'] = cv_results['estimator'][0].feature_names_in_
            df_run['Score'] = cv_run.feature_importances_
            df_run['CV_run'] = [index+1] * len(cv_run.feature_importances_)
            df_run['Average_score'] = feature_importances
            self.selected_features_scores = pd.concat([self.selected_features_scores, df_run], ignore_index=True)

        # drop features that were not selected and sort by average importance score
        self.selected_features_scores = self.selected_features_scores[self.selected_features_scores['Feature'].isin(self.selected_features_names)]
        self.selected_features_scores.sort_values('Average_score', ascending=False, inplace=True)

        # print results
        print(f"\tCV test scores feature selection (R2): {cv_results['test_score']}")

        # save selected features in model
        self.X_train_selected = self.X_train[self.X_train.columns.intersection(self.selected_features_names)]
        self.X_test_selected = self.X_test[self.X_test.columns.intersection(self.selected_features_names)]

    def select_features_by_genetic_algorithm(self):
        """
        Select feature by genetic algorithm using genetic_selection.GeneticSelectionCV
        todo: visualize scores by generation (if available in gsv._support), accept parameters
        todo: @jose not sure if these default parameters are ok, what do you think?
        """
        print(f'-> select features with genetic algorithm')
        # Create GeneticSelectionCV to optimize feature selection using genetic algorithm
        gscv = GeneticSelectionCV(self.regressor, cv=5, scoring='r2', n_population=100,
                                  n_gen_no_change=5, n_jobs=-1, max_features=30, caching=True, verbose=True)
        gscv.fit(self.X_train, self.y_train.values)

        # # Get the selected features
        selected_features_bool = gscv.get_support()
        self.selected_features_names = self.features_names[selected_features_bool]
        self.X_train_selected = self.X_train.iloc[:, selected_features_bool]
        self.X_test_selected = self.X_test[self.X_test.columns.intersection(self.selected_features_names)]

        self.settings_string += "_GeneticSelection"

        # save cv importance scores for visualization
        for cv_run in [1,2,3,4,5]: # 5 inner CV runs
            df_run = pd.DataFrame()
            df_run['Feature'] = self.selected_features_names
            df_run['Score'] = [np.nan] * len(self.selected_features_names)
            df_run['CV_run'] = [cv_run] * len(self.selected_features_names)
            df_run['Average_score'] = [np.nan] * len(self.selected_features_names)
            self.selected_features_scores = pd.concat([self.selected_features_scores, df_run], ignore_index=True)

    def select_features_by_sequential_selector(self, method ='top_k', top_k = 50, min_improvement = 0.001, verbose = False):
        """
        Select features by using SequentialFeatureSelection in forward direction.

        @param method: 'top_k' (default) or 'min_improvement'
        @param top_k: top k features to select
        @param min_improvement: minimal improvement of the score when adding a feature
        """
        # check input validity and save parameters for later
        print(f'-> Select features by Sequential Feature Selector.')

        if method == 'top_k':
            # first, check if enough features available for feature selection by top_k
            if not self.has_enough_features_for_selection(top_k):
                return

            print(f'-> Selecting top {top_k} features')
            self.settings_string += f"_k{top_k}"
            n_features = len(self.X_train.columns)
            if top_k > n_features:
                top_k = n_features
                print('Warning: top_k = {} features to be selected, but only {} features available:'.format(top_k, n_features))
                print('--> top_k was set to {} '.format(n_features))
        elif method == 'min_improvement':
            print(f'-> Selecting features that improve the score by more than {min_improvement} (max features: top_k = {top_k})')
            self.settings_string += f"_min{min_improvement}"
        else:
            raise ValueError("The 'method' parameter should be set to 'top_k' or 'min_improvement'")

        # run forward SFS
        sfs = SequentialFeatureSelector(self.regressor, cv=5, scoring='r2', n_jobs=-1, k_features=top_k)
        sfs.fit(self.X_train, self.y_train.values)

        # save cv importance scores for visualization
        current_average_score = current_min_score = current_max_score = -1000 # set to a really low number to ensure that at least one feature is selected
        previous_set = set([])
        for i in range(1,top_k+1): # loop through increments to fetch features and scores
            df_run = pd.DataFrame()
            iteration = sfs.subsets_[i]
            new_set = set(iteration['feature_idx'])
            new_feature_index = new_set.difference(previous_set).pop()
            previous_set = new_set

            # check if there is no more improvement in average, min and max cv score by adding a new feature
            improvement =  ((np.mean(iteration['cv_scores']) - current_average_score) > min_improvement) \
                or ((np.min(iteration['cv_scores']) - current_min_score) > min_improvement) \
                or ((np.max(iteration['cv_scores']) - current_max_score) > min_improvement)
            if verbose:
                print(f"CV scores: {iteration['cv_scores']}, improvement={improvement}")

            if method == 'min_improvement' and not improvement:
                break
            else:
                # save new scores
                current_average_score = np.mean(iteration['cv_scores'])
                current_min_score = np.min(iteration['cv_scores'])
                current_max_score = np.max(iteration['cv_scores'])
                feature_name = self.features_names[new_feature_index]

                df_run['Feature'] = [feature_name] * sfs.cv
                df_run['Score'] = iteration['cv_scores']
                df_run['CV_run'] = range(1,sfs.cv+1)
                df_run['Average_score'] = [current_average_score] * sfs.cv
                self.selected_features_scores = pd.concat([self.selected_features_scores, df_run], ignore_index=True)
                self.selected_features_names.append(feature_name)

        # Get the selected features
        self.X_train_selected = self.X_train[self.X_train.columns.intersection(self.selected_features_names)]
        self.X_test_selected = self.X_test[self.X_test.columns.intersection(self.selected_features_names)]

    def has_enough_features_for_selection(self, top_k):
        """
        Checks if enough features are available in X_train for feature selection (top_k)
        @param top_k: number of features to be selected
        @return: False if more features to be selected than available, True otherwise
        """
        n_features = len(self.X_train.columns)
        if top_k >= n_features:
            print('Warning: top_k = {} features to be selected, but only {} features available:'.format(top_k,
                                                                                                        n_features),
                  '--> Feature selection is skipped.')
            # set selected features to all features
            self.no_feature_selection()
            self.feature_selection_method = None # set feature selection to None
            return False
        else:
            return True

    def reduce_feature_space(self):
        """
        Dimensionality reduction of feature space. Supports PCA and SVD
        """
        kwargs = self.feature_reduction_params
        if self.feature_reduction_method == 'pca':
            self.reduce_feature_dimensionality(method_name='PCA', **kwargs)
        elif self.feature_reduction_method == 'svd':
            self.reduce_feature_dimensionality(method_name='SVD', **kwargs)
        elif self.feature_reduction_method in [None, 'None']:
            self.no_feature_reduction()
        else:
            raise ValueError(
                "Possible values for feature_reduction methods: 'pca', 'svd', None")

    def no_feature_reduction(self):
        """
        If no feature dimensionality reduction is applied,
        X_train/test_reduced are directly set to X_train/test_selected
        """
        print(f'-> No feature dimensionality reduction applied')
        self.X_train_reduced = self.X_train_selected
        self.X_test_reduced = self.X_test_selected
        self.reduced_features_names = self.X_train_selected.columns

    def reduce_feature_dimensionality(self, n_components, method_name):
        """
        Reduce feature space using Principal Component Analysis (PCA)
        :param n_components: number of components to be used as features, or list vales\
        :param method_name: PCA or SVD
        """
        print(f'-> Reduce feature space with {method_name} to {n_components} components')
        n_features = len(self.X_train.columns)
        if n_components > n_features:
            n_components = n_features
            print('Warning: n_components = {} to be generated, but only {} features available:'.format(n_components,
                                                                                                        n_features))
            print('--> n_components was set to {} '.format(n_features))

        if method_name == 'PCA':
            method = PCA()
        elif method_name == 'SVD':
            method = TruncatedSVD()
        else:
            raise ValueError(f"Method {method_name} should be either 'PCA' or 'SVD'")

        # if several values for n_components are given: find optimal number of components
        if type(n_components) == list:
            if len(n_components) != 1:
                # perform grid search
                pipe = Pipeline([("reduce_dim", "passthrough"), ("regressor", self.regressor)] )
                param_grid = [
                    {
                        "reduce_dim": [method],
                        "reduce_dim__n_components": n_components
                    }
                ]
                grid = GridSearchCV(pipe, cv=5, n_jobs=1, param_grid=param_grid, return_train_score=True)
                # get best components
                grid.fit(self.X_train_selected, self.y_train.values)
                n_components = grid.best_params_['reduce_dim__n_components']
                print("\tOptimal number of components:", n_components)
                print("\tAverage scores:", grid.cv_results_['mean_test_score'])
            else:
                n_components = n_components[0]

        self.settings_string += f"_PC{n_components}"
        method.set_params(n_components=n_components)

        # perform PCA on all X_train_selected
        method.fit(self.X_train_selected, self.y_train.values)
        self.X_train_reduced = method.transform(self.X_train_selected)
        if not self.X_test_selected.empty:
            self.X_test_reduced = method.transform(self.X_test_selected)
        self.feature_transformer = method

        print(method.explained_variance_ratio_)
        self.reduced_features_names = []
        [self.reduced_features_names.append(f'PC{i+1}') for i in range(0,n_components)]


    ############## Splitting #########################

    def basic_random_split(self, state, test_size=None):
        """Splits the data randomly and assigns the necessary attributes for training and testing.
        Note that to optimize the models only the training set is used, and it is further split using cross validation.
        The test set is never seen during optimization.
        The concept is similar to nested cross validation except that this first split is random.
        For a typical nested_cross_validation see define_outer_loop()

        Parameters
        ----------
        state: The random state used to split the data.
        test_size: Default is 0.2 (i.e., 20%)
        """
        if test_size is None:
            test_size = self.test_size

        self.X_train, self.X_test, self.y_train, self.y_test, self.y_train_std, self.y_test_std = train_test_split(
                                                                                self.features,
                                                                                self.target_variable,
                                                                                self.target_variable_std,
                                                                                test_size=test_size,
                                                                                random_state=state)
        # transform y
        self.y_train = self.y_train[self.target_variable_name]
        self.y_test = self.y_test[self.target_variable_name]
        self.test_indices = self.y_test.index.values.tolist()

    def define_outer_loop(self, train_index, test_index):
        """Splits the data using the  provided train and test indices.
        This method is expected to be used as part of nested_cross_val_screening()
        which does a nested cross validation for several (regressor, features) pairs.
        The idea is to screen which pairs seem to perform better on a given dataset.

        Parameters
        ----------
        train_index: The indices for the training set.
        test_index: The indices for the test set.
        """
        self.train_indices = train_index
        self.test_indices = test_index
        self.X_train = self.features.iloc[train_index]
        self.X_test = self.features.iloc[test_index]
        self.y_train = self.target_variable.iloc[train_index][self.target_variable_name]
        self.y_test = self.target_variable.iloc[test_index][self.target_variable_name]
        if self.has_target_variable_std():
            self.y_train_std = self.target_variable_std.iloc[train_index]
            self.y_test_std = self.target_variable_std.iloc[test_index]

    def use_all_data_for_training(self):
        """
        Use all data for training - no test set is created
        """
        train_indices = self.features.index.tolist()
        test_indices = []
        self.define_outer_loop(train_indices, test_indices)

    def split_by_compound(self, state=0, test_size=None):
        """Uses the static method 'random_smiles_split' to split the data "by_compound" meaning that it will create
        a list of unique SMILES strings then split that list according to the desired test_size.
        The lists are used to create dataframes with the required information to later get X_train, y_train, etc.

        This function is necessary when the input data includes the same compound multiple times
        (e.g., using values of the same compound but in different treatment plants or different conditions)
        but you want to make sure that the model is tested using only previously unseen molecules.

        Parameters
        ----------
        state: The random state used to split the data.
        test_size: Default is 0.2 (i.e., 20%)
        """
        if test_size is None:
            test_size = self.test_size

        smiles_train, smiles_test = self.random_smiles_split(list(self.data[self.smiles_name].unique()),
                                                             split_size=test_size, random_state=state)

        train_data = self.data[self.data[self.smiles_name].isin(smiles_train)]
        test_data = self.data[self.data[self.smiles_name].isin(smiles_test)]

        # The idea is to select the training set based on smiles_name but later avoid SMILES as part of the input
        train_subset = self.features[self.features[self.smiles_name].isin(smiles_train)]
        self.X_train = train_subset.loc[:, train_subset.columns != self.smiles_name]

        # The idea is to select the training set based on smiles_name but later avoid SMILES as part of the input
        test_subset = self.features[self.features[self.smiles_name].isin(smiles_test)]
        self.X_test = test_subset.loc[:, test_subset.columns != self.smiles_name]

        self.y_train = train_data[self.target_variable_name]
        self.y_test = test_data[self.target_variable_name]
        self.test_indices = test_data.index.values.tolist()
        return


    @staticmethod
    def random_smiles_split(smiles_list, split_size=0.2, random_state=None):
        """
        Receives a list (SMILES in this case) as input, shuffles it, splits it according to split_size, and returns
        two lists (for train and test sets).
        This function is used to make sure that compounds appear only once in train or test sets.

        :param smiles_list: List of unique smiles strings; typically from the curated data.
        :param split_size: Fraction intended to be allocated to the test set
        :param random_state:
        :return:
        """
        if random_state is not None:
            random.seed(random_state)

        # Calculate the size of each split
        split_point = int(len(smiles_list) * split_size)

        # Randomly shuffle the input list
        shuffled_list = random.sample(smiles_list, len(smiles_list))

        # Split the shuffled list into two parts
        train_smiles = shuffled_list[split_point:]
        test_smiles = shuffled_list[:split_point]

        return train_smiles, test_smiles

    def regression_stratified_split(self):
        binsplits = range(0, 25, 3)
        y_categorized = pd.cut(self.target_variable, bins=binsplits,
                               include_lowest=True, right=False, labels=range(0, 24, 3))
        print(y_categorized)
        # todo:  Ignore for now. The idea is to create equally spaced "categories"
        #  of a continuous variable and then sample from each category
        #  so that the train and test sets contain examples along the whole range of target_variable_values

    ######################### Training and testing #########################################

    def train_and_test(self, verbose=True):
        """
        Simple function to train and test a regression model.
        @param verbose: print output
        """
        if verbose:
            print("\n############# Training ############# ")

        self.regressor.fit(self.X_train, self.y_train)

        self.y_pred_train = self.regressor.predict(self.X_train)
        self.y_pred_train = pd.DataFrame(self.y_pred_train)
        self.y_pred_train.index = self.X_train.index
        if verbose:
            print("\n############# Testing ############# ")
        self.y_pred = self.regressor.predict(self.X_test)
        self.y_pred = pd.DataFrame(self.y_pred)
        self.y_pred.index = self.X_test.index

    def simple_run(self,
                   random_state_list=None,
                   verbose=True,
                   split_by_compound=False,
                   test_size=None):
        """
        Models are created and validated based on a random split of available data.
        If a random_state_list is provided it will rerun and validate for each random split.
        Training and testing scores are printed.
        Train score refers to metrics based on "predicting" with the training set.
        The purpose is to understand how well the model fits the training data.
        Test score refers to metrics based on actual predictions on previously unseen data.
        These scores are stored as attributes of model as "test_scores_df".
        Note that it prints average statistics

        :param verbose: Show random state
        :param random_state_list:
        :param split_by_compound: If True it makes sure that compounds appear only once in the training or test sets.
        It is necessary when showing the same substances more than once (e.g., same substance in different conditions)
        but you want to test on unseen molecules.
        :param test_size: if None it takes the default 20%.
        """
        if random_state_list is None:
            random_state_list = [0]
        for random_state in random_state_list:
            if verbose:
                print("Random state: {}".format(random_state))

            if split_by_compound:
                self.split_by_compound(state=random_state)
                print('split was done by compound')

            else:
                self.basic_random_split(random_state, test_size=test_size)
                print('split was done completely random')

            self.train_and_test(verbose=verbose)
            self.save_fitting_values(random_state)
            self.save_predicted_values(random_state)
            self.save_scores(random_state, 'train', verbose=verbose)
            self.save_scores(random_state, 'test', verbose=verbose)

        # Print average scores
        self.print_average_statistics()

    def simple_cross_val_run(self, run_id, verbose=True):
        """
        Models are created and validated based on a random split of available data.
        If a random_state_list is provided it will rerun and validate for each random split.
        Training and testing scores are printed.
        Train score refers to metrics based on "predicting" with the training set.
        The purpose is to understand how well the model fits the training data.
        Test score refers to metrics based on actual predictions on previously unseen data.
        These scores are stored as attributes of model as "test_scores_df".
        :param verbose: Show random state
        :param run_id: The fold number to identify the run; it is important when analyzing the scores.
        """

        if verbose:
            print("Random state: {}".format(random_state))

        self.train_and_test(verbose=verbose)
        self.save_fitting_values(run_id)
        self.save_predicted_values(run_id)
        self.save_scores(run_id, 'train', verbose=verbose)
        self.save_scores(run_id, 'test', verbose=verbose)

    def train_regressor_on_subset(self, subset_index: []):
        """
        Train the regressor on a subset of the data, defined by subset_index
        @param subset_index: indices of data subset for model training
        """
        print(f'-> Train {self.regressor_name} regressor on subset')

        # reduce training
        self.X_train = self.X_train.loc[subset_index]
        self.y_train = self.y_train.loc[subset_index]
        if self.has_target_variable_std():
            self.y_train_std = self.y_train_std.loc[subset_index]

        # select features from preprocessed and pre-selected X_train
        self.select_features() # select features based on regressor
        self.reduce_feature_space() # apply dimensionality reduction

        # train final model()
        self.set_regressor_parameters(self.regressor, self.regressor_params)
        self.train_model()


    def size_dependent_run(self, random_state_list=None, verbose=True):
        """
        Models are created and validated based on a random split of available data.
        If a random_state_list is provided it will rerun and validate for each random split.
        Training and testing scores are printed.
        Train score refers to metrics based on "predicting" with the training set.
        The purpose is to understand how well the model fits the training data.
        Test score refers to metrics based on actual predictions on previously unseen data.
        These scores are stored as attributes of model as "test_scores_df".
        :param verbose: Show random state
        :param random_state_list:
        """
        if random_state_list is None:
            random_state_list = [0]
        for random_state in random_state_list:

            # Perform initial split into validation and test sets
            X_remaining, self.X_test, y_remaining, self.y_test = train_test_split(self.features,
                                                                                  self.target_variable,
                                                                                  test_size=self.test_size,
                                                                                  random_state=random_state)
            self.test_indices = self.y_test.index.values.tolist()
            # Dictionary to store results
            results = {}

            # Loop through different fractions of validation set
            for train_fraction in np.arange(0.1, 1.1, 0.1):
                # Calculate the size of training set based on the fraction of the remaining data
                train_index = int(train_fraction * len(X_remaining))

                # Extract the training data from the remaining data
                self.X_train, self.y_train = X_remaining[:train_index], y_remaining[:train_index]

                self.train_and_test(verbose=verbose)
                self.save_fitting_values(random_state)
                self.save_predicted_values(random_state)
                self.save_scores(random_state, 'train', verbose=verbose, train_fraction=train_fraction)
                self.save_scores(random_state, 'test', verbose=verbose, train_fraction=train_fraction)

        # Print average scores
        self.print_average_statistics()
        return

    ######################### Saving methods ##############################################

    def save_fitting_values(self, random_state):
        """ Must be called after train to get fitted values"""
        self.fitting_data_values = pd.DataFrame()
        self.fitting_data_values['modeled'] = self.y_pred_train
        self.fitting_data_values['experimental'] = self.y_train
        self.fitting_data_values['random_state'] = random_state
        # self.fitting_data_values.to_csv(self.fitting_data_tsv, sep='\t', index=True)
        self.training_dict[random_state] = self.fitting_data_values

    def save_predicted_values(self, random_state):
        """ Must be called after train and test to get predicted values"""
        self.predicted_target_variable = pd.DataFrame()
        self.predicted_target_variable['SMILES'] = self.data[self.smiles_name][self.test_indices]
        self.predicted_target_variable['predicted'] = self.y_pred
        self.predicted_target_variable['experimental'] = self.y_test
        self.predicted_target_variable['random_state'] = random_state
        self.predicted_target_variable['absolute_error'] = (
            abs(self.predicted_target_variable['predicted'] - self.predicted_target_variable['experimental']))

        if self.has_y_pred_score:
            self.predicted_target_variable['predicted_score'] = self.y_pred_score
        if self.has_target_variable_std():
            self.predicted_target_variable['experimental_std'] = self.y_test_std

        output_filename = os.path.join(self.get_data_directory(),
                                 'predictions_{}_{}_{}_{}.csv'.format(self.data_type, self.tag,
                                                                              self.setup_name, self.settings_string))
        print("-> save predictions to", output_filename)
        self.predicted_target_variable.to_csv(output_filename)

    def save_external_predictions(self, descriptors):
        """ Same as save_predicted_values, but for external predictions. """
        self.predicted_target_variable = pd.DataFrame()
        self.predicted_target_variable[self.smiles_name] = descriptors.features[self.smiles_name]
        self.predicted_target_variable[self.target_variable_name + '_predicted'] = self.y_pred
        if self.has_y_pred_score:
            self.predicted_target_variable[self.target_variable_std_name + '_predicted'] = self.y_pred_score

        # add experimental data if available from self.data
        self.predicted_target_variable = pd.merge(self.predicted_target_variable,
                                                  self.data.loc[:, self.data.columns != self.id_name],
                                                  on=self.smiles_name, how='left')
        # rename columns to distinguish experimental from predicted values
        self.predicted_target_variable.rename(columns={self.target_variable_name: self.target_variable_name+'_experimental',
                                                       self.target_variable_std_name: self.target_variable_std_name+'_experimental'},
                                              inplace=True)


    def save_scores(self, run_id: int = 0, stage: str = 'test', train_fraction: float = 1, verbose=True):
        """
        Save training and testing scores.

        @param run_id: id to identify for which split the results were obtained.
        It can be a random state of a fold number during kfold
        @param stage: 'test' or 'train'
        @param train_fraction: Fraction of training data used for this specific run
        @param verbose: verbose
        @return:
        """
        if stage == 'train':
            y_true = self.y_train
            y_pred = self.y_pred_train

        elif stage == 'test':
            y_true = self.y_test
            y_pred = self.y_pred

        else:
            y_true = None
            y_pred = None

        scores_dic = {
            'r2': [r2_score(y_true, y_pred)],
            'mse': [mean_squared_error(y_true, y_pred, squared=True)],
            'rmse': [mean_squared_error(y_true, y_pred, squared=False)],
            'mae': [mean_absolute_error(y_true, y_pred)],
            'descriptors': [self.descriptors.get_current_feature_space()],
            'regressor': [self.regressor_name],
            'feature_selection': [str(self.feature_selection_method)],
            'feature_reduction': [str(self.feature_reduction_method)],
            'run_id': [run_id],
            'train_fraction': [train_fraction],
        }

        if stage == 'train':
            self.train_scores = pd.DataFrame(scores_dic)
            if verbose:
                print("Train scores: \n{}".format(self.train_scores.to_markdown()))
            self.train_scores_df = pd.concat([self.train_scores_df, self.train_scores], ignore_index=True)
        elif stage == 'test':
            self.test_scores = pd.DataFrame(scores_dic)
            if verbose:
                print("Test scores: \n{}\n{}".format(self.test_scores.to_markdown(), '-'*50))
            self.test_scores_df = pd.concat([self.test_scores_df, self.test_scores], ignore_index=True)
        else:
            print("Stage not clearly defined")
            return None

    def print_average_statistics(self):
        print("\n############# Statistics ############# ")
        print("Regressor: {}\nFeature space: {}".format(str(self.regressor), self.descriptors.current_feature_space))
        print('Target variable length: {}\nNumber of features: {}\nUncertainty metric: {}'.format(
            self.target_variable.size,
            self.features.shape[1],
            self.target_variable_std_name))
        print("Average train scores:")
        print(self.train_scores_df.loc[:, ['r2', 'mse', 'rmse', 'mae']].mean())
        print("Average test scores:")
        print(self.test_scores_df.loc[:, ['r2', 'mse', 'rmse', 'mae']].mean())



    def perform_grid_search(self, X, save_results = False):
        print("-> optimize hyperparameters:", self.regressor_params)
        function_name = sys._getframe().f_code.co_name  # get the name of the function
        # ensure that self.regressor_params has multiple options

        # define pipeline
        pipe = Pipeline([
            ('regressor', self.regressor)  # Regressor
        ])
        # pipeline parameters need to start with 'regressor__'
        pipeline_params = {}
        for p in self.regressor_params:
            pipeline_params['regressor__' + p] = self.regressor_params[p]

        # Create GridSearchCV to optimize hyperparameters
        grid_search = GridSearchCV(pipe, pipeline_params, cv=5, scoring='r2')
        grid_search.fit(X, self.y_train.values)
        print("\tCV test scores hyperparameter tuning (R2):", grid_search.cv_results_['mean_test_score'])

        # save optimized parameters
        self.best_regressor_model = grid_search.best_estimator_['regressor']
        self.best_regressor_params = grid_search.best_estimator_['regressor'].get_params()

        # save to file
        if save_results:
            self.save_optimized_settings(function_name)

        print(f"Final hyperparameters: {self.best_regressor_params}")

    def predict_target_variable(self, descriptors: Descriptor):
        """
        Predict target variable for a descriptors object
        @param descriptors: Descriptor object for input SMILES
        """
        # retrieve only features that are used by the model
        X = descriptors.features[descriptors.features.columns.intersection(self.feature_names_used_for_training)]

        # perform min-max scaling
        print("-> preprocess features (min-max scaling)")
        temp_feature_names = self.feature_preprocessor.feature_names_in_
        feature_dict = {}
        for name in temp_feature_names:
            if name in X.columns:
                feature_dict[name] = X[name]
            else:
                feature_dict[name] = 0
        temp_features = pd.DataFrame(feature_dict)
        X_preprocessed = self.feature_preprocessor.transform(temp_features)

        # if a transformer for feature dimensionality reduction is available
        if self.feature_transformer:
            X_selected = X_preprocessed[X_preprocessed.columns.intersection(self.feature_transformer.feature_names_in_)]
            X_reduced = self.feature_transformer.transform(X_selected)
        else: # if no dimensionality reduction was performed
            X_selected = X_preprocessed[X_preprocessed.columns.intersection(self.feature_names_used_for_training)]
            X_reduced = X_selected

        # run prediction
        self.predict(X_reduced)
        # save output
        self.save_external_predictions(descriptors)

    def predict(self, X):
        if self.regressor_name == 'Gaussian Process Regressor':
            self.y_pred, self.y_pred_score = self.regressor.predict(X, return_std=True)
        elif self.regressor_name == "KNN Regressor":
            self.y_pred = self.regressor.predict(X)
            # get average distance to closest neighbors
            neighbors = self.regressor.kneighbors(X,n_neighbors=self.regressor_params['n_neighbors'], return_distance=True)
            scores = neighbors[0][:,0] # get only distance of the nearest neighbor
            # scores = np.mean(neighbors[0], axis=1) # other option: get average scores over all neighbors,
            # normalize scores (distances) between 0 and 1
            self.y_pred_score = (scores-np.min(scores))/(np.max(scores)-np.min(scores))
        else:
            self.y_pred = self.regressor.predict(X)


    def complete_train_regressors(self):
        """

        @param feature_selection: Type of feature selection 'importance', 'genetic', 'sequential'
        """
        print(f"Regressor: {self.regressor_name}")
        function_name = sys._getframe().f_code.co_name  # get the name of the function

        # optimize hyperparameters
        self.perform_grid_search(self.X_train, save_results=True)

        # Select features
        self.select_features()  # Here features are truly selected self.features_names -> self.selected_features_names
        self.reduce_feature_space() # Here the space is reduced

        # optimize hyperparameters after feature selection, but only if features were selected or reduced
        if self.feature_selection_method not in [None, 'None'] or self.feature_reduction_method not in [None, 'None']:
            self.perform_grid_search(self.X_train_reduced, save_results=True)

        # config and train final model()
        self.set_regressor_parameters(self.best_regressor_model, self.best_regressor_params)
        self.train_model()

        # visualize selected features where applicable
        if self.feature_selection_method in ['sequential', 'importance']:
            v = Visualize(self, function_name)
            v.feature_selection_plot()

    def set_regressor_parameters(self, regressor, parameters: dict):
        """
        Set fixed regressor parameters when no optimization is performed
        @param regressor: Regressor object
        @param parameters: parameter dictionary
        """
        self.regressor = regressor
        self.regressor_params = {}
        for p in parameters:
            if type(parameters[p]) == list:
                assert len(parameters[p]) == 1, f'Several parameters are given for {p} but no grid search is performed'
                self.regressor_params[p] = parameters[p][0]
            else:
                self.regressor_params[p] = parameters[p]
        self.regressor.set_params(**self.regressor_params)

    def train_model(self):
        """
        Fit self.regressor with self.X_train_reduced and self.y_train
        """
        if self.regressor_name == 'Gaussian Process Regressor' and self.has_target_variable_std():
            self.regressor.set_params(alpha=self.y_train_std.values)
        self.regressor.fit(self.X_train_reduced, self.y_train.values)

    def save_optimized_settings(self, analysis_type):
        """
        Save optimized settings (regressor, parameters, features) to .txt file
        @param analysis_type: name of analysis type used for output file name
        """
        # save settings
        self.build_optimized_settings_files(analysis_type=analysis_type)
        open_file = open(self.optimized_settings_txt, 'w')
        open_file.write('Best regressor model: {}\n'.format(self.best_regressor_model))
        open_file.write('Best regressor parameters: {}\n'.format(self.best_regressor_params))
        open_file.write('Selected features names: {}\n'.format(self.selected_features_names))
        open_file.close()

        # save file with selected features scores
        self.selected_features_scores.to_csv(self.selected_features_scores_tsv, sep='\t')


    def build_optimized_settings_files(self, analysis_type):
        """
        Build filenames to save optimized settings as .txt files
        @param analysis_type: name of analysis type used for output file name
        """
        file_dir = os.path.join(self.get_data_directory(), 'optimized_settings', analysis_type, self.curation_type)
        self.build_directory_structure(file_dir)
        if not self.optimized_settings_txt:
            self.optimized_settings_txt = self.build_output_filename(
                os.path.join(file_dir, 'optimized_settings_{}'.format(self.regressor_name)), suffix='.txt')
        if not self.selected_features_scores_tsv:
            self.selected_features_scores_tsv = self.build_output_filename(
                os.path.join(file_dir, 'selected_features_scores{}'.format(self.regressor_name)))


    def get_train_data_fractions(self, fraction_list):
        """
        Obtain a dictionary of indices of training data for each fraction in fraction_list
        @param fraction_list: list of fractions of training data to be used for training, e.g., [0.1, 0.2, ..., 1]
        @return: dictionary of indices of training data, where keys are fractions from input fraction_list
        """
        shuffled_indices = self.train_indices
        random.shuffle(shuffled_indices)
        fraction_dict = {}
        for fraction in fraction_list:
            stop_index = int(fraction * len(shuffled_indices))
            # Extract the training data from the remaining data
            indices = shuffled_indices[:stop_index]
            indices.sort()
            fraction_dict[np.round(fraction,1)] = indices
        return fraction_dict

    def complete_evaluate_models(self, run_id, train_fraction = 1, visualize=True):
        """
        Evaluate a model by predicting y for the test set, saving scores, and visualizing predicted vs. true values
        @param run_id: split ID, e.g. from k-fold CV
        @param train_fraction: fraction of training data that was used for training. By default, all training data is used
        @param visualize:
        """
        function_name = sys._getframe().f_code.co_name

        # predict target values for test set
        self.predict_for_test_set()

        # save scores to files
        self.save_scores(run_id, 'test', train_fraction)
        self.save_predicted_values(run_id)

        # visualize prediction on test set
        if visualize:
            v = Visualize(self, function_name)
            v.scatterplot_predicted_vs_test()
            # Additional plots for regressors providing prediction scores
            if self.regressor_name == 'Gaussian Process Regressor':
                threshold_list = [0.5, 0.65, 0.8, 1] # for GPR, this represents the predicted standard deviation of y
                v.plot_performance_vs_score_threshold(threshold_list)
                v.scatterplots_by_thresholds(threshold_list)
            elif self.regressor_name == 'KNN Regressor':
                threshold_list = [0.2, 0.4, 0.6, 1] # for KNN, this is the average distance to neighbors
                v.plot_performance_vs_score_threshold(threshold_list) # , reverse = True
                v.scatterplots_by_thresholds(threshold_list)  # , reverse = True

    def predict_for_test_set(self):
        """
        Predict the target variable (y_pred) and, if available, a prediction score (y_pred_score) for the test set
        """
        # predict y for test set
        self.predict(self.X_test_reduced)

    def drop_missing_target_variable_value(self):
        """
        The purpose of this method is to drop compounds with missing data for the target variable.
        As an example, in the case of WWTP data, I want to calculate descriptors for all the relevant compounds;
        yet, in some cases there is influent and/or effluent data missing for specific treatment plants.
        """

        print("\n############# Drop compounds without data for this specific plant ############# ")
        print("Shape before dropping missing endpoint is:{}".format(self.data.shape))
        self.data.replace(r'', np.nan, regex=True)
        # self.data.replace(-np.inf, np.nan, regex=True)  # todo: check that this works with snf_dat
        self.data.dropna(axis='index', subset=self.target_variable_name, inplace=True)
        print("Shape after dropping missing endpoint is:{}".format(self.data.shape))

    def update_modeling_data(self):
        """
        If samples are removed from model data AFTER the descriptors have been calculated
        then self.features and self.model_data will not match in length so the purpose of this method
        is to merge them based on index and redefine the features and target variables.
        """
        joint_data = pd.merge(self.features, self.data, on=self.id_name)

        self.data = joint_data.loc[:, [self.id_name,
                                       self.smiles_name,
                                       self.compound_name,
                                       self.target_variable_name]]

        self.features = joint_data[self.features.columns]
        self.target_variable = self.data[self.target_variable_name]
        return
