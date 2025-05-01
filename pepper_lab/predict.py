import os
import joblib

import numpy as np
import pandas as pd
from rdkit import Chem

from pepper_lab.descriptors import Descriptors
from pepper_lab.pepper import Pepper
from pepper_lab.util import Util


class Predict(Pepper):
    def __init__(self, renku=False):
        """
        Initiate Predict object
        :param renku: set to True if the predictions run on renku
        """
        super().__init__()
        pep = Pepper(renku=renku)
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

        :param input_model: Pepper.Model object
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
            self.model = Predict.load_pickle(input_model)
        self.tag = self.model.tag
        self.data_type = self.model.data_type

        # load smiles
        self.set_smiles_name(self.model.smiles_name)
        self.check_smiles_input(input_smiles, input_smiles_type)

        # calculate descriptors and predict endpoints
        if not self.descriptors.model_data.empty: # if at least some of the smiles are valid
            self.descriptors.set_smiles_name(self.model.smiles_name)
            self.descriptors.set_data_type(self.data_type) # get the data type from the model used for prediction
            print(self.descriptors.smiles_name)
            self.descriptors.load_descriptors(from_csv=precalculated_descriptors, load_by_feature_name=True,
                                              feature_name_list = self.model.feature_names_used_for_training,
                                              feature_space_map=self.model.descriptors.feature_space_map)
            # fetch feature space from original model and set it for descriptors of external data
            self.descriptors.define_feature_space(self.model.descriptors.get_current_feature_space())
            self.model.predict_target_variable(self.descriptors, use_individual_trees=self.model.use_individual_trees)
            output_df = self.create_prediction_output_table()
        else: # if no valid smiles
            output_df = self.input_data
            output_df[self.target_variable_name + '_predicted'] = np.nan
            output_df[self.target_variable_std_name + '_predicted'] = np.nan

        #save to file
        output_file_path = self.build_output_filename('output/'+self.file_tag)
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
        return output_df

    def check_smiles_input(self, input_smiles, input_smiles_type):
        print("-> checking SMILES input")
        if input_smiles_type == 'tsv':
            path_to_file = self.data_directory + '/input/' + input_smiles
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
    def load_pickle(input_model):
        return joblib.load(open(input_model, 'rb'))
