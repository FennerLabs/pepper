"""
This example script reproduces our initial attempts with default regressors.
The differences with the original example script (predict_breakthrough_wwtp.py) are the following:

- No additional curation of training data
- default hyperparameters
- Simple 5-fold cross-validation (as opposed to nested cross-validation)
- complete chemical space (['padel', 'mordred', 'maccs', 'mfps', 'rdkitfps', 'ep_trig'])
- Load descriptors from_csv=True (added them manually to avoid fetching ep_trig from envipath which takes time)
- All regressors (['MLR', 'Ridge', 'KRidge', 'SGD', 'LSVR', 'RF', 'GB', 'AB', 'DT', 'MLP', 'SVR', 'KNN'])
"""

from pepper_lab.pepper import Pepper
from pepper_lab.datastructurewwtp import DataStructureWWTP
from pepper_lab.descriptors import Descriptors
from pepper_lab.modeling import Modeling
from pepper_lab.visualize import Visualize


if __name__ == '__main__':
    # ------------------------------------#
    # Initialize with user preferences  #
    # ------------------------------------#
    pep = Pepper(renku=True)
    pep.set_setup_name('example')
    pep.set_data_type('WWTP')
    pep.set_tag('combined_data')
    pep.set_target_variable_name('logB')
    pep.set_target_variable_std_name('stdev')
    pep.set_smiles_name('CanonicalSMILES')

    # ------------------------------------#
    # Data import, curation and analysis #
    # ------------------------------------#
    wwtp_data = DataStructureWWTP(pep)
    print(wwtp_data.name_data_tsv)
    wwtp_data.load_data('name_data')
    wwtp_data.load_data('full_data')
    wwtp_data.set_curation_type('basic_curation_nitrification_plants_only_newtest20250701')
    wwtp_data.set_plant_list(wwtp_data.ndn_plants)
    wwtp_data.curate_by_target_variable()
    wwtp_data.set_id_name('Combined_ID')
    wwtp_data.reduce_for_modelling(no_curation=True)
    wwtp_data.select_modeling_data()
    wwtp_data.create_modelling_input()
    print("Data loaded and curated.")

    # ------------------------------------#
    # Descriptor calculation              #
    # ------------------------------------#
    # calculate descriptors
    descriptors = Descriptors(pep)
    descriptors.set_data(wwtp_data)
    descriptors.load_descriptors(from_csv=True,
                                 PaDEL=True, mordred=True,
                                 MACCS=True, enviPath_trig=True,
                                 mfps=True, RDKit_fps=True)
    
    print("Descriptors loaded.")

    # ------------------------------------#
    # Modeling                            #
    # ------------------------------------#

    wwtp_modeling = Modeling(pep, wwtp_data, descriptors)
    print("Modeling initialized.")

    # wwtp_modeling.set_regressor_settings(mode='singlevalue', config='default_wwtp')
    wwtp_modeling.test_different_setups_cross_val(regressor_name_list = ['MLR', 'Ridge', 'KRidge', 'SGD', 'LSVR', 'RF', 'GB', 'AB', 'DT', 'MLP', 'SVR', 'KNN'],
                                            feature_space_list = ['padel', 'mordred', 'maccs', 'mfps', 'rdkitfps', 'ep_trig'],
                                            # feature_space_list = ['maccs'],
                                            mode='singlevalue',
                                            config='default_wwtp')
    
    print('done')
    