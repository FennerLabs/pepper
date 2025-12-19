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
# set false for reproducing results from scratch and build final model
reproduce_results = False

# Set true to run TP prediction at the end
TP_prediction = True

# Set to True to run zeroPM prediction at the end
zeropm_prediciton = False



if __name__ == "__main__":
    print(os.path.abspath(os.sep))
    current_path = os.getcwd()
    pep = Pepper(pepper_data_location=os.path.join(current_path, "../.."))
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
    soil_descriptors.load_descriptors(from_csv=not reproduce_results, enviPath_prob=False, enviPath_trig=True, mordred=False,
                                       PaDEL=True, MACCS=True, RDKit_descriptors=True, avalonfps=True, RDKit_fps=True)
#    ----------------- 
#    ---Modeling------
#    -----------------
    soil_data.analyze_target_variable_distributions()
    soil_modeling = Modeling(pep, soil_data, soil_descriptors)
    
    if reproduce_results:
        soil_modeling.nested_cross_validation_screening(regressor_name_list = ['GPR', 'RF'],
                                                feature_space_list=['ep_trig','padel', "maccs", 'rdkitfps', 'avalonfps', 'rdkitdesc','all'], load_existing=False, config= 'soil_paper') 


    # build final model with best hyperparemeters for GPR

    soil_model = soil_modeling.build_final_model(regressor_name='GPR', feature_space='padel', config = 'soil_paper_GPR_optimized')


    if TP_prediction:
        ############
        #Prediciton#
        ############

        new_prediction = Predict(pep)

        # load TPs without DT50 (downloaded from www.enviPath.com; accessed: 12.04.25)
        TP_without = pd.read_csv("../data/soil/TP_soil_without_DT50.tsv", sep='\t')
        # load training data
        parent_compounds = pd.read_csv("../data/soil/full_data_soil_all_data.tsv", sep='\t')

        # use the same SMILES processing as in datastructure.py
        TP_without['smiles'] = TP_without['smiles'].apply(lambda x: Util.canonicalize_smiles(
                        Util.remove_stereo_info(Util.remove_stereochemistry(x))))

        # filter out TPs without DT50 that also are part of the training data
        parent_set = set(parent_compounds[['compound_name', 'SMILES']].apply(tuple, axis=1))
        mask = (
        TP_without['compound_name'].isin(parent_compounds['compound_name'])
        | TP_without['smiles'].isin(parent_compounds['SMILES']))
        filtered = TP_without[~mask]
        filtered.rename(columns={'smiles': 'SMILES'}, inplace=True)
        filtered.drop_duplicates(subset=['SMILES'], inplace=True)

        # create SMILES input file for prediction
        filtered['SMILES'].to_csv(os.path.join(new_prediction.get_data_directory(), "input/TP_without_DT50_smiles.tsv"), sep='\t')
        
        # calculate descriptors for TP smiles and predict
        input_file_name = "TP_without_DT50_smiles.tsv"

        # use the pickled model
        pickled_model_filename = os.path.join(soil_modeling.get_data_directory(),'models/' f'final_model_{soil_model.regressor_name_short}.pkl')
        new_prediction.predict_endpoint(input_model=pickled_model_filename,input_model_format='pickle', input_smiles=input_file_name, input_smiles_type='tsv', precalculated_descriptors=not reproduce_results)

        ###----------------------------------------------
        ###---Visualization of TP and parent compounds---
        ###----------------------------------------------

        ### load data ###

        final_predictions = new_prediction.model.predicted_target_variable.copy()
        final_predictions.to_csv(new_prediction.data_directory + '/output/' + 'final_predictions_soil_without_TP_raw.tsv', sep='\t', index=False)
        print(f"final predictions: average logDT50_mean:{final_predictions['logDT50_mean_predicted'].mean()}\n average logDT50_std:{final_predictions['logDT50_std_predicted'].mean()}\n number of compounds: {final_predictions.shape[0]}" )

        # calculate probabilities for final predictions
        final_predictions['p_nP'], final_predictions['p_P'], final_predictions['p_vP'] = zip(
        *final_predictions.apply(lambda row: new_prediction.model.calculate_prediction_probabilities(row['logDT50_mean_predicted'], row['logDT50_std_predicted'], T_P = 120, T_vP = 180), axis=1))
        final_predictions['p_nP'] = (final_predictions['p_nP'] * 100).round(0).astype(int)
        final_predictions['p_P'] = (final_predictions['p_P'] * 100).round(0).astype(int)
        final_predictions['p_vP'] = (final_predictions['p_vP'] * 100).round(0).astype(int)
        final_predictions['logDT50_mean_predicted'] = final_predictions['logDT50_mean_predicted'].round(2)
        final_predictions['logDT50_std_predicted'] = final_predictions['logDT50_std_predicted'].round(2)

        # load TP meta data
        TP_meta_data = pd.read_csv("../data/soil/TP_soil_without_DT50.tsv", sep='\t')
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
        parent_compounds = pd.read_csv(new_prediction.data_directory + '/../data_structure/soil/model_data_soil_all_data.tsv', sep='\t')
        # load raw data
        raw_data = pd.read_csv(os.path.join(new_prediction.data_directory, '../data_structure/soil/cpd_data_soil_all_data.tsv'), sep='\t')
        
        # filter for all TPs that are in training data
        training_TPs = raw_data[raw_data['node_depth'] != 0]

        # save training TPs
        training_TPs.to_csv(new_prediction.data_directory + '/output/' + 'training_TP_soil_all_data.tsv', sep='\t', index=False)


        # merge TP meta data pathway information with final predictions
        overlapping_smiles = set(final_predictions['SMILES']).intersection(set(TP_meta_data['SMILES']))
        final_predictions = final_predictions[final_predictions['SMILES'].isin(overlapping_smiles)]
        final_predictions = final_predictions.merge(TP_meta_data[['SMILES', 'pathway','minormajor', 'compound_name']], on='SMILES', how='left')
        final_predictions.rename(columns={'compound_name_y': 'compound_name'}, inplace=True)
        final_predictions.drop(columns=['compound_name_x'], inplace=True)

        # save final predictions file
        final_predictions.to_csv(new_prediction.data_directory + '/output/' + 'Predictions_pesticide_TPs.csv', index=False)
        
        
        #TODO: finish for reproducing plots
        
        """
        # drop duplicates
        # TP_meta_data.drop_duplicates(subset=['compound_name'], inplace=True)
        raw_data.drop_duplicates(subset=['compound_name', 'pathway_name'], inplace=True)

        # filter our parent compounds
        depth_information = raw_data[raw_data['node_depth'] == 0]

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
        save_path_pc = new_prediction.data_directory + '/output/' + 'parent_compounds_with_class_probabilities.tsv'
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
        save_path_TP = new_prediction.data_directory + '/output/' + 'Predictions_pesticide_TPs_filtered.tsv'
        final_predictions_plot.to_csv(save_path_TP,sep="\t", index=False) 

        final_predictions_plot['minormajor'] = final_predictions_plot['minormajor'].fillna('')

        print(len(final_predictions))
        final_predictions_plot.drop(columns=['SMILES', 'minormajor']).to_latex(new_prediction.data_directory + '/output/' + 'final_predictions_with_persistence_probs_filterd.tex', index=False, float_format="%.2f")

        """
    if zeropm_prediciton:

        ###------------------------
        ###---Zero PM prediction---
        ###------------------------

        # Set of inorganic elements and metals
        INORGANIC_ATOMS = set([
            "Na","K","Li","Mg","Ca","Fe","Cu","Zn",
            "Al","Mn","Co","Ni","Cr","Pb","Hg"
        ])


        def is_organic(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
            
            # must contain carbon
            if "C" not in atoms:
                return False
            
            # no metals/inorganics
            if any(a in INORGANIC_ATOMS for a in atoms):
                return False
            
            return True
        
        def is_not_salt(smiles):
            return "." not in smiles
        
        def mw_below_1200(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            mw = dp.MolWt(mol)
            return mw < 1200
        
        def smiles_passes_filters(smiles):
            return (
                is_not_salt(smiles) and
                is_organic(smiles) and
                mw_below_1200(smiles)
            )
        
        input_file_name = "raw_zeropm.csv"
        zeropm = pd.read_csv(os.path.join(pep.get_data_directory(), 'predict/', 'input/', input_file_name))

        
        print("Number of zeroPM:",len(zeropm))
        num_salts = sum(zeropm["SMILES"].apply(lambda s: "." in s))
        print("Number of salts:", num_salts)
        num_over_1200 = sum(zeropm["SMILES"].apply(lambda s: not mw_below_1200(s)))
        print("Number of molecules > 1200 Da:", num_over_1200)
        filtered_df = zeropm[zeropm["SMILES"].apply(smiles_passes_filters)]

        filtered_df['SMILES'].to_csv(os.path.join(pep.get_data_directory(), "predict/input/filtered_zeropm_smiles.tsv"), sep='\t', index=False)

        # #TODO: do not run this since padel descriptors are for rest zero pm only
        # zeropm_file_csv= pd.read_csv(os.path.join(pep.get_data_directory(), "predict/input/filtered_zeropm.tsv"), sep='\t')
        # zeropm_file_csv['SMILES'].to_csv(os.path.join(pep.get_data_directory(), "predict/input/filtered_zeropm_smiles.tsv"), sep='\t', index=False)
        
        # zero_pm_prediction = Predict()
        # zero_pm_prediction.data_type = 'soil_zeropm'
        # zeropm_file = "filtered_zeropm_smiles.tsv"  
        # zero_pm_prediction.predict_endpoint(input_model=pickled_model_filename,input_model_format='pickle', input_smiles=zeropm_file, input_smiles_type='tsv', precalculated_descriptors=True)

        # zero_pm_predictions = zero_pm_prediction.model.predicted_target_variable.copy()



        # print(f"Zero PM predictions: average logDT50_mean:{zero_pm_predictions['logDT50_mean_predicted'].mean()}\n average logDT50_std:{zero_pm_predictions['logDT50_std_predicted'].mean()}")
        # print(f"final predictions: average logDT50_mean:{final_predictions['logDT50_mean_predicted'].mean()}\n average logDT50_std:{final_predictions['logDT50_std_predicted'].mean()}")

        # # filter out parent compounds
