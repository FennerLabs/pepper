"""
------------
---Import---
------------
"""
# python 
import sys
import os
import random

# third-party
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns   
import enviPath_python as ep

# project
from pepper_lab.pepper import *
from pepper_lab.datastructure import *
from pepper_lab.descriptors import *
from pepper_lab.datastructuresoil import *
from pepper_lab.modeling import *
from pepper_lab.predict import *

"""
----------
---Main---
----------
"""

# set to True for reproducing results from scratch and build final model
reproduce_results = False

# Set to True to run predictions for TPs
TP_prediction = False

# Set to True to run predictions for zeroPM substances
zeropm_prediction = False

# Set to True to run predictions for benchmarking compounds
benchmarking_prediction = False


if __name__ == "__main__":
    print(os.path.abspath(os.sep))
    current_path = os.getcwd()
    pep = Pepper(pepper_data_location=os.path.join(current_path, "..", ".."))
    pep.set_tag('all_data')
    pep.set_setup_name('uncertainty_paper_soil')
    pep.set_data_type('soil')
    pep.set_target_variable_name('logDT50_mean')
    pep.set_target_variable_std_name('logDT50_std')
    pep.set_smiles_name('SMILES')
    pep.set_compound_name('compound_name')

   
#   -----------------
#   ---Data Import---
#    -----------------
    soil_data = DataStructureSoil(pep)
    soil_data.curate_annotate(from_csv=True, from_paper=True)
    soil_data.reduce_for_modelling(from_csv=False)

#    -----------------
#    ---Descriptors---
#    -----------------
    soil_descriptors = Descriptors(pep)
    soil_descriptors.set_data(soil_data)
    soil_descriptors.load_descriptors(from_csv=False, enviPath_prob=False, enviPath_trig=True, mordred=False,
                                       PaDEL=True, MACCS=True, RDKit_descriptors=True, avalonfps=True, RDKit_fps=True)
#    ----------------- 
#    ---Modeling------
#    -----------------
#     soil_data.analyze_target_variable_distributions()
    soil_modeling = Modeling(pep, soil_data, soil_descriptors)
    
    if reproduce_results:
        soil_modeling.nested_cross_validation_screening(regressor_name_list = ['GPR', 'RF'],
                                                        feature_space_list=['ep_trig','padel', "maccs", 'rdkitfps',
                                                                            'avalonfps', 'rdkitdesc','all'],
                                                        load_existing=False, config= 'soil_paper')

    # build final model with best hyperparameters for GPR
    soil_model = soil_modeling.build_final_model(regressor_name='GPR', feature_space='padel',
                                                 config = 'soil_paper_GPR_optimized')

    # save model file name for further use
    pickled_model_filename = os.path.join(soil_modeling.get_data_directory(),
                                          'models', f'final_model_{soil_model.regressor_name_short}.pkl')


    if TP_prediction:
        ############
        #Prediction#
        ############

        new_prediction = Predict(pep)

        # load TPs without DT50 (downloaded from www.enviPath.com; accessed: 12.04.25)
        path_to_TPs_without_DT50 = os.path.join("..", "data", "soil", "TP_soil_without_DT50.tsv")
        TP_without = pd.read_csv(path_to_TPs_without_DT50, sep='\t')
        # load training data
        training_compounds = soil_data.full_data

        # use the same SMILES processing as in datastructure.py
        TP_without['smiles'] = TP_without['smiles'].apply(lambda x: Util.canonicalize_smiles(
                        Util.remove_stereo_info(Util.remove_stereochemistry(x))))

        # filter out TPs without DT50 that also are part of the training data
        parent_set = set(training_compounds[['compound_name', 'SMILES']].apply(tuple, axis=1))
        mask = (
                TP_without['compound_name'].isin(training_compounds['compound_name'])
                | TP_without['smiles'].isin(training_compounds['SMILES'])
        )
        filtered = TP_without[~mask]
        filtered.rename(columns={'smiles': 'SMILES'}, inplace=True)
        filtered.drop_duplicates(subset=['SMILES'], inplace=True)

        # create SMILES input file for prediction and save to input folder
        input_file_name = "TP_without_DT50_smiles.tsv"
        filtered['SMILES'].to_csv(os.path.join(new_prediction.get_data_directory(), "input", input_file_name), sep='\t')

        # calculate descriptors for TP smiles and predict endpoint
        new_prediction.predict_endpoint(input_model=pickled_model_filename,input_model_format='pickle',
                                        input_smiles=input_file_name, input_smiles_type='tsv',
                                        precalculated_descriptors=False) # set to true if descriptors are already calculated

        ###----------------------------------------------
        ###---Visualization of TP and parent compounds---
        ###----------------------------------------------

        ### load data ###

        final_predictions = new_prediction.model.predicted_target_variable.copy()
        final_predictions.to_csv(os.path.join(new_prediction.data_directory,
                                              'output',
                                              'final_predictions_soil_without_TP_raw.tsv'),
                                 sep='\t', index=False)
        print(f"final predictions: "
              f"\n\taverage logDT50_mean: {round(final_predictions['logDT50_mean_predicted'].mean(), 2)}"
              f"\n\taverage logDT50_std: {round(final_predictions['logDT50_std_predicted'].mean(), 2)}"
              f"\n\tnumber of compounds: {final_predictions.shape[0]}" )

        # calculate probabilities for final predictions
        final_predictions['p_nP'], final_predictions['p_P'], final_predictions['p_vP'] = zip(
        *final_predictions.apply(
            lambda row:
            new_prediction.model.calculate_prediction_probabilities(row['logDT50_mean_predicted'],
                                                                    row['logDT50_std_predicted'],
                                                                    T_P = 120, T_vP = 180), axis=1))
        final_predictions['p_nP'] = (final_predictions['p_nP'] * 100).round(0).astype(int)
        final_predictions['p_P'] = (final_predictions['p_P'] * 100).round(0).astype(int)
        final_predictions['p_vP'] = (final_predictions['p_vP'] * 100).round(0).astype(int)
        final_predictions['logDT50_mean_predicted'] = final_predictions['logDT50_mean_predicted'].round(2)
        final_predictions['logDT50_std_predicted'] = final_predictions['logDT50_std_predicted'].round(2)

        # load TP meta data
        TP_meta_data = pd.read_csv(path_to_TPs_without_DT50, sep='\t')
        TP_meta_data.rename(columns={'smiles': 'SMILES'}, inplace=True)

        # filter out co2 smiles
        TP_meta_data = TP_meta_data[TP_meta_data['compound_name'].str.contains('CO2') == False]
        final_predictions = final_predictions[final_predictions[new_prediction.smiles_name] != 'O=C=O']

        # filter out water and other single ions
        exclude_list = [
        "[Zn-]", "[SH-]", "[O-][Cl+][O-]", "[NH4+]", 
        "[Cl-]", "[CH3+]", "[Br-]", "O"]

        final_predictions = final_predictions[~final_predictions["SMILES"].isin(exclude_list)]

        # load parent compound data
        parent_compounds = soil_data.model_data
        # load raw data
        full_data = soil_data.full_data
        
        # filter for all TPs that are in training data
        training_TPs = full_data[full_data['node_depth'] != 0]

        # save training TPs
        training_TPs.to_csv(os.path.join(new_prediction.data_directory, 'output', 'training_TP_soil_all_data.tsv'),
                            sep='\t', index=False)


        # merge TP meta data pathway information with final predictions
        overlapping_smiles = set(final_predictions['SMILES']).intersection(set(TP_meta_data['SMILES']))
        final_predictions = final_predictions[final_predictions['SMILES'].isin(overlapping_smiles)]
        final_predictions = final_predictions.merge(TP_meta_data[['SMILES', 'pathway','minormajor', 'compound_name']],
                                                    on='SMILES', how='left')
        final_predictions.rename(columns={'compound_name_y': 'compound_name'}, inplace=True)
        final_predictions.drop(columns=['compound_name_x'], inplace=True)

        # save final predictions file
        final_predictions.to_csv(os.path.join(new_prediction.data_directory, 'output', 'Predictions_pesticide_TPs.csv'), index=False)

        # drop duplicates
        # TP_meta_data.drop_duplicates(subset=['compound_name'], inplace=True)
        full_data.drop_duplicates(subset=['compound_name', 'pathway_name'], inplace=True)

        # filter out parent compounds
        depth_information = full_data[full_data['node_depth'] == 0]

        # merge parent compounds with raw data to get pathway information
        parent_compounds = soil_data.model_data.merge(full_data[['pathway_name' , 'compound_name']],
                                            left_on='compound_name', right_on='compound_name', how='left')
        
        # merge parent compounds with depth information to get node depth
        parent_compounds = parent_compounds.merge(depth_information[['compound_name', 'node_depth']],
                                            left_on='compound_name', right_on='compound_name', how='left')

        pipu = parent_compounds[parent_compounds['node_depth'] == 0]['pathway_name'].value_counts()

        # get rid of the ones with more than 1 compound from the parent compounds and keep only the ones with one compound
        parent_compounds = parent_compounds[parent_compounds['pathway_name'].isin(pipu[pipu == 1].index)]

        # sort parent compounds by logdt50 mean

        # parent_compounds['logDT50_mean'] = parent_compounds.groupby('pathway_name')['logDT50_mean'].transform('max')
        parent_compounds.sort_values(by='logDT50_mean', ascending=True, inplace=True)

        # Filter final predictions for std_prediction 
        print("min and max predicted logDT50_mean:", final_predictions['logDT50_mean_predicted'].min(), final_predictions['logDT50_mean_predicted'].max())
        print("min and max predicted logDT50_std:", final_predictions['logDT50_std_predicted'].min(), final_predictions['logDT50_std_predicted'].max())

        final_predictions = final_predictions[final_predictions['logDT50_std_predicted'] < 0.7]

        # Keep only pathways that are in both parent_compounds and final_predictions
        valid_pathways = parent_compounds['pathway_name'].unique()
        final_predictions = final_predictions[final_predictions['pathway'].isin(valid_pathways)]

        # Filter parent_compounds accordingly
        parent_compounds = parent_compounds[parent_compounds['pathway_name'].isin(final_predictions['pathway'].unique())]

        # Define consistent pathway order

        # Use only parent compounds (node_depth == 0) to define ordering
        parent_only = parent_compounds[parent_compounds['node_depth'] == 0].copy()

        # Sort parents by their logDT50_mean
        parent_only.sort_values(by='logDT50_mean', ascending=True, inplace=True)

        # Define ordered pathways based on parent compounds only
        ordered_pathways = parent_only['pathway_name'].tolist()

        # Apply categorical ordering
        final_predictions['pathway'] = pd.Categorical(final_predictions['pathway'], categories=ordered_pathways, ordered=True)
        parent_compounds['pathway_name'] = pd.Categorical(parent_compounds['pathway_name'], categories=ordered_pathways, ordered=True)

        # filter out parent compounds pathways that are not in final_predictions
        parent_compounds = parent_compounds[parent_compounds['pathway_name'].isin(final_predictions['pathway'].unique())]

        v = Visualize(new_prediction.model, analysis_type='soil')

        v.plot_predicted_compounds(final_predictions=final_predictions, parent_compounds=parent_compounds, ordered_pathways=ordered_pathways)

        
        ####-----------------
        # Calculate class probabilities
        ####-----------------

        # filter first s.t table matches plot
        parent_compounds_plot = parent_compounds[(parent_compounds['node_depth'] == 0) &
        (parent_compounds['pathway_name'].isin(final_predictions['pathway'].unique()))
        ].copy()

        # calculate probabilities for filtered parent compounds
        parent_compounds_plot['p_nP'], parent_compounds_plot['p_P'], parent_compounds_plot['p_vP'] = zip(
        *parent_compounds_plot.apply(lambda row: new_prediction.model.calculate_prediction_probabilities(row['logDT50_mean'], row['logDT50_std'], T_P = 120, T_vP = 180), axis=1)
        )

        # ensure ordering matches plot
        parent_compounds_plot['pathway_name'] = pd.Categorical(
        parent_compounds_plot['pathway_name'],
        categories=ordered_pathways,
        ordered=True
        )

        parent_compounds_plot['logDT50_mean'] = parent_compounds_plot['logDT50_mean'].round(2)
        parent_compounds_plot['logDT50_std'] = parent_compounds_plot['logDT50_std'].round(2)

        parent_compounds_plot['p_nP'] = (parent_compounds_plot['p_nP'] * 100).round(0).astype(int)
        parent_compounds_plot['p_P'] = (parent_compounds_plot['p_P'] * 100).round(0).astype(int)
        parent_compounds_plot['p_vP'] = (parent_compounds_plot['p_vP'] * 100).round(0).astype(int)  

        # save to tsv
        save_path_pc = os.path.join(new_prediction.data_directory, 'output', 'parent_compounds_with_class_probabilities.tsv')
        parent_compounds_plot.to_csv(save_path_pc, sep='\t', index=False)

        # same for final predictions

        final_predictions_plot = final_predictions.copy()

        # Match order to ordered_pathways
        final_predictions_plot['pathway'] = pd.Categorical(
        final_predictions_plot['pathway'],
        categories=ordered_pathways,
        ordered=True
        )

        # Sort by pathway category order and save to tsv
        final_predictions_plot = final_predictions_plot.drop(columns=['logDT50_mean_experimental', 'logDT50_std_experimental'])
        final_predictions_plot = final_predictions_plot.sort_values('pathway', kind='stable')
        save_path_TP = os.path.join(new_prediction.data_directory, 'output', 'Predictions_pesticide_TPs_filtered.tsv')
        final_predictions_plot.to_csv(save_path_TP,sep="\t", index=False) 

        final_predictions_plot['minormajor'] = final_predictions_plot['minormajor'].fillna('')

        print("Number of TPs with predictions", len(final_predictions))
        # Optional: Save as latex table. Requires Jinja2 to be installed
        # figure_predictions = os.path.join(new_prediction.data_directory, 'output', 'final_predictions_with_persistence_probs_filterd.tex')
        # final_predictions_plot.drop(columns=['SMILES', 'minormajor']).to_latex(figure_predictions, index=False, float_format="%.2f")


    if zeropm_prediction:

        ###------------------------
        ###---Zero PM prediction---
        ###------------------------

        def is_organic(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
            # must contain carbon
            if "C" not in atoms:
                return False
            return True
        
        def is_not_salt(smiles):
            return "." not in smiles
        
        def mw_below_1200(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            mw = Chem.Descriptors.MolWt(mol)
            return mw < 1200
        
        def smiles_passes_filters(smiles):
            return (
                is_not_salt(smiles) and
                is_organic(smiles) and
                mw_below_1200(smiles)
            )

        # Data: ZeroPM database of marketed chemicals
        # The data was downloaded via https://pubchem.ncbi.nlm.nih.gov/source/25168 on 10/10/2025
        # Query for ZeroPM substances :
        # --> https://pubchem.ncbi.nlm.nih.gov//#query=Sg3v6Dv9XkFpa1xy3goVVMyEgOSUWaHQ2_W6nMDkqJ3A_ZQ&selected_id_type=cid
        # Query for ZeroPM compounds :
        # --> https://www.ncbi.nlm.nih.gov/pccompound?cmd=HistorySearch&hinit=true&query_key=8&WebEnv=MCID_68e8e425746347fdcf94f5e6
        # The 140,161 ZeroPM substances were linked to 137,915 unique PubChem compounds (due to duplicate CID entries in ZeroPM list)
        # The downloaded file PubChem_compound_cache_46RGTSiBTT16E8UKR3KMLFXwgJC3R6Z-3Fu9MsdKrzPHU5M.csv
        # was renamed to raw_zeropm.csv and saved under pepper_data/predict/input/

        input_file_name = "raw_zeropm.csv"
        # input_file_name = "raw_zeropm_subset.csv" # use a subset to test the code
        zeropm = pd.read_csv(os.path.join(pep.get_data_directory(), 'predict', 'input', input_file_name))

        print("Number of zeroPM:",len(zeropm))
        num_salts = sum(zeropm["SMILES"].apply(lambda s: "." in s))
        print("Number of salts:", num_salts)
        num_over_1200 = sum(zeropm["SMILES"].apply(lambda s: not mw_below_1200(s)))
        print("Number of molecules > 1200 Da:", num_over_1200)
        filtered_df = zeropm[zeropm["SMILES"].apply(smiles_passes_filters)]

        filtered_df['SMILES'].to_csv(os.path.join(pep.get_data_directory(), "predict", "input", "filtered_zeropm_smiles.tsv"), sep='\t', index=False)

        zero_pm_prediction = Predict(pep)
        zero_pm_prediction.data_type = 'soil_zeropm'
        zeropm_file = "filtered_zeropm_smiles.tsv"
        zero_pm_prediction.predict_endpoint(input_model=pickled_model_filename, input_model_format='pickle',
                                            input_smiles=zeropm_file, input_smiles_type='tsv',
                                            precalculated_descriptors=False) # set to true if already calculated

        zero_pm_predictions = zero_pm_prediction.model.predicted_target_variable.copy()

        print(f"Zero PM predictions: "
              f"\n\taverage logDT50_mean: {round(zero_pm_predictions['logDT50_mean_predicted'].mean(), 2)}"
              f"\n\taverage logDT50_std: {round(zero_pm_predictions['logDT50_std_predicted'].mean(), 2)}")


    if benchmarking_prediction:
        # setup pepper object
        pep_benchmark = Pepper(pepper_data_location=os.path.join(current_path, "..", ".."))
        pep_benchmark.set_tag('25') # 25 substances
        pep_benchmark.set_setup_name('benchmarking')
        pep_benchmark.set_data_type('soil')
        pep_benchmark.set_target_variable_name('logDT50_mean')
        pep_benchmark.set_target_variable_std_name('logDT50_std')
        pep_benchmark.set_smiles_name('SMILES')
        pep_benchmark.set_compound_name('compound_name')

        # create benchmarking data structure
        benchmarking_data = DataStructureSoil(pep_benchmark)
        benchmarking_data.load_data('raw_data', source='data')
        benchmarking_data.curate_annotate(from_csv=False, curate=False)
        benchmarking_data.reduce_for_modelling(from_csv=False, curate=False)

        # predict half-lives
        new_prediction = Predict(pep)
        new_prediction.set_data_type('soil')
        new_prediction.set_tag('benchmarking')
        new_prediction.predict_endpoint(input_model = soil_model,
                                        input_smiles=benchmarking_data.model_data,
                                        input_smiles_type='dataframe')
