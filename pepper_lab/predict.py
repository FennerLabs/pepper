import os
import joblib

import numpy as np
import pandas as pd
from rdkit import Chem

from pepper_lab.descriptors import Descriptors
from pepper_lab.pepper import Pepper
from pepper_lab.util import Util
from pepper_lab.visualize import Visualize


class Predict(Pepper):
    def __init__(self, pep: Pepper, renku=False):
        """
        Initiate Predict object
        :param renku: set to True if the predictions run on renku
        """
        super().__init__()
        pep = pep
        pep.set_tag('prediction')
        self.set_data_directory(os.path.join(pep.data_directory, 'predict'))
        self.build_directory_structure('input')
        self.build_directory_structure('output')
        self.descriptors = Descriptors(pep)
        self.input_data = pd.DataFrame()
        self.file_tag = ''  # used to designate input file (for .tsv input) and output file
        self.model = None  # model used for predictions


    def predict_endpoint(self, input_model, input_smiles, input_model_format='model',
                         input_smiles_type='tsv', precalculated_descriptors = False):
        """
        Given a Pepper.Model object and a (list of) SMILES, this function predicts endpoints for the input structure(s)
        and saves the prediction (incl. experimental data points, if available) to the pepper_data/predict/output/ folder.

        :param input_model: Pepper.Model object or path to a pickle file containing the model object
        :param input_smiles: input SMILES as str objects, pandas DataFrame, or input file name (tab_seperated, to be
        saved under pepper_data/predict/inputt. Mandatory columnn is 'SMILES', other columns are optional and will be
        copied to the output file.
        :param input_model_format: 'model' or 'pickle'
        :param input_smiles_type: 'tsv' (tab-separated txt file) or 'smi' (e.g., 'c1ccccc1') or 'dataframe' (column header must match pepper.smiles_name)
        :param precalculated_descriptors: set to true when descriptors are provided
        """
        print('\n############# Predict endpoints ############# ')
        # load model
        if input_model_format == 'model':
            self.model = input_model
        elif input_model_format == 'pickle':
            self.model = Predict.load_joblib(input_model)
        self.tag = self.model.tag
        self.data_type = ""
        #self.data_type = self.model.data_type
        self.model.prediction_mode = True  # set model to prediction mode

        # load smiles
        self.set_smiles_name(self.model.smiles_name)
        self.check_smiles_input(input_smiles, input_smiles_type)

        valid_input = True
        # calculate descriptors and predict endpoints
        if self.descriptors.model_data.empty:
            valid_input = False
        else: # if at least some of the smiles are valid
            self.descriptors.set_smiles_name(self.model.smiles_name)
            self.descriptors.set_data_type(self.data_type) # get the data type from the model used for prediction
            self.descriptors.load_descriptors(from_csv=precalculated_descriptors, load_by_feature_name=True,
                                              feature_name_list = self.model.feature_names_used_for_training,
                                              feature_space_map=self.model.descriptors.feature_space_map)
            # fetch feature space from original model and set it for descriptors of external data
            try:
                self.descriptors.define_feature_space(self.model.descriptors.get_current_feature_space())
            except AssertionError: # feature space was empty
                valid_input = False

        if valid_input: # claculate model predictions
            self.model.predict_target_variable(self.descriptors, use_individual_trees=self.model.use_individual_trees)
            self.model.create_prediction_probabilities()
            output_df = self.create_prediction_output_table()

        else: # if no valid smiles or features, return empty df
            output_df = self.input_data
            output_df[self.target_variable_name + '_predicted'] = np.nan
            output_df[self.target_variable_std_name + '_predicted'] = np.nan

        #save to file
        output_file_path = self.build_output_filename(os.path.join("output", 'predicted'))
        output_df.to_csv(output_file_path, sep='\t', index=False)
        print("Predictions are saved to {}".format(output_file_path))

        return output_df

    def create_prediction_output_table(self): # if predictions only
        """
        Create output table with comments
        """
        print('-> create prediction output table')
        # Add comments
        value_counts = self.input_data[self.smiles_name].value_counts()
        new_warning_list = []
        # collect warnings
        for index, row in self.input_data.iterrows():
            warning = row['warnings']
            if row[self.smiles_name] in self.model.predicted_target_variable[self.smiles_name].values:
                if value_counts.get(row[self.smiles_name] , 0) > 1:
                    warning += 'compound duplicated in input file'
            else:
                if warning != '':
                    warning += ', '
                warning += 'descriptors could not be calculated'
            new_warning_list.append(warning)
        self.input_data['warnings'] = new_warning_list

        output_df = self.input_data.merge(self.model.predicted_target_variable,
                                          on=self.smiles_name, how='left')  # get predictions + scores where we have them.
        
        output_df = output_df.merge(self.model.prediction_probabilities, on=self.smiles_name, how='left')  # get class probabilities where we have them.

        # # remove columns that are not needed
        # output_df.drop(columns=['logDT50_mean_experimental','logDT50_std_experimental','compound_name'], inplace=True, errors='ignore')
        
        return output_df

    def check_smiles_input(self, input_smiles, input_smiles_type):
        print("-> checking SMILES input")
        if input_smiles_type == 'tsv':
            path_to_file = os.path.join(self.data_directory, 'input', input_smiles)
            self.descriptors.model_data = pd.read_csv(path_to_file, sep='\t', encoding_errors='ignore')
            self.file_tag = input_smiles.split('.')[0]
        elif input_smiles_type == 'smi':
            self.descriptors.model_data = pd.DataFrame({self.smiles_name: [input_smiles]})
            self.file_tag = 'single_smiles'
        elif input_smiles_type == 'dataframe':
            self.descriptors.model_data = input_smiles
            self.file_tag = 'dataframe'
        else:
            raise NotImplementedError("Please provide valid input smiles type")

        checked_smiles = []
        warnings = []
        for smiles in self.descriptors.model_data[self.smiles_name]:
            try:
                mol = Chem.MolFromSmiles(smiles, sanitize=True)
            except Exception as e:
                print('Text: {} \n not recognized as a SMILES string'.format(e))
                mol = None

            if mol is None:
                warnings.append('SMILES not valid')
                checked_smiles.append(np.nan)
                print('SMILES not valid:', smiles)
            else:
                can = Util.canonicalize_smiles(smiles)
                checked_smiles.append(can)
                warnings.append('')

        # Data to keep
        self.input_data['original_' + self.smiles_name] = self.descriptors.model_data[self.smiles_name]
        for column_name in [self.id_name, self.compound_name]:
            if column_name in self.descriptors.model_data.columns:
                self.input_data[column_name] = self.descriptors.model_data[column_name]

        # new data generated
        self.input_data[self.smiles_name] = checked_smiles
        self.input_data['warnings'] = warnings
        self.descriptors.model_data[self.smiles_name] = checked_smiles
        self.descriptors.model_data.dropna(axis='rows', inplace=True)

    @staticmethod
    def load_joblib(input_model):
        return joblib.load(open(input_model, 'rb'))
    

    def calculate_class_probabilities(self):
        pass

    def visualize_predictions(self):
        final_predictions = self.model.predicted_target_variable.copy()

 
        ### load data ###
        # load TP meta data
        TP_meta_data = pd.read_csv(self.data_directory + '/input/' + 'TP_soil_without_DT50.tsv', sep='\t')
        TP_meta_data.rename(columns={'smiles': 'SMILES'}, inplace=True)
        # load parent compound data
        parent_compounds = pd.read_csv(self.data_directory + '/../data_structure/soil/model_data_soil_all_data.tsv', sep='\t')
        # load raw data
        raw_data = pd.read_csv(os.path.join(self.data_directory, '../data_structure/soil/raw_data_soil_all_data.tsv'), sep='\t')
        # load depth information
        depth_information = pd.read_csv(self.data_directory + '/input/' + 'cpd_data_soil_all_data_parent_node_depth.tsv', sep='\t')

        # merge TP meta data pathway information with final predictions
        TP_meta_data['SMILES'] = TP_meta_data['SMILES'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x), isomericSmiles=True))
        overlapping_smiles = set(final_predictions['SMILES']).intersection(set(TP_meta_data['SMILES']))
        final_predictions = final_predictions[final_predictions['SMILES'].isin(overlapping_smiles)]
        final_predictions = final_predictions.merge(TP_meta_data[['SMILES', 'pathway']], on='SMILES', how='left')
      
        # filter out co2 smiles
        TP_meta_data = TP_meta_data[TP_meta_data['compound_name'].str.contains('CO2') == False]
        final_predictions = final_predictions[final_predictions[self.smiles_name] != 'O=C=O']

        # drop duplicates
        TP_meta_data.drop_duplicates(subset=['compound_name'], inplace=True)
        raw_data.drop_duplicates(subset=['compound_name', 'pathway_name'], inplace=True)

        # filter our parent compounds
        depth_information = depth_information[depth_information['node_depth'] == 0]

        # merge parent compounds with raw data to get pathway information
        parent_compounds = parent_compounds.merge(raw_data[['pathway_name' , 'compound_name']],
                                          left_on='compound_name', right_on='compound_name', how='left')
        
        # merge parent compounds with depth information to get node depth
        parent_compounds = parent_compounds.merge(depth_information[['compound_name', 'node_depth']],
                                          left_on='compound_name', right_on='compound_name', how='left')

        pipu = parent_compounds[parent_compounds['node_depth'] == 0]['pathway_name'].value_counts() 

        # get rid of the ones with more than 1 compound from the parent compounds and keep only the ones with one compound
        parent_compounds = parent_compounds[parent_compounds['pathway_name'].isin(pipu[pipu == 1].index)]

        # sort parent compounds by logdt50 mean
        parent_compounds['logDT50_mean'] = parent_compounds.groupby('pathway_name')['logDT50_mean'].transform('max')
        parent_compounds.sort_values(by='logDT50_mean', ascending=True, inplace=True)

        # Filter final predictions for std_prediction < 0.6
        final_predictions = final_predictions[final_predictions['logDT50_std_predicted'] < 0.6]

        # Keep only pathways that are in both parent_compounds and final_predictions
        valid_pathways = parent_compounds['pathway_name'].unique()
        final_predictions = final_predictions[final_predictions['pathway'].isin(valid_pathways)]

        # Filter parent_compounds accordingly
        parent_compounds = parent_compounds[parent_compounds['pathway_name'].isin(final_predictions['pathway'].unique())]

        # Define consistent pathway order
        ordered_pathways = parent_compounds['pathway_name'].drop_duplicates().tolist()

        # Apply categorical ordering
        final_predictions['pathway'] = pd.Categorical(final_predictions['pathway'], categories=ordered_pathways, ordered=True)
        parent_compounds['pathway_name'] = pd.Categorical(parent_compounds['pathway_name'], categories=ordered_pathways, ordered=True)

        # filter out parent compounds pathways that are not in final_predictions
        parent_compounds = parent_compounds[parent_compounds['pathway_name'].isin(final_predictions['pathway'].unique())]

        v = Visualize(self.model, analysis_type='soil')

        v.plot_predicted_compounds(final_predictions=final_predictions, parent_compounds=parent_compounds)