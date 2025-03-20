# This is an example ML workflow using pepper-lab

from pepper.pepper import Pepper
# from pepper-lab.datastructure import DataStructure
from pepper.datastructurewwtp import DataStructureWWTP
from pepper.descriptors import Descriptors
from pepper.modeling import Modeling


if __name__ == '__main__':
    # ------------------------------------#
    # Initialize with user preferences  #
    # ------------------------------------#
    pep = Pepper()  # Set renku=True in case of running a Renku session
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
    wwtp_data.load_data('name_data', source='data')
    wwtp_data.load_data('full_data', source='data')
    wwtp_data.set_curation_type('testing_pypi_version')
    wwtp_data.set_plant_list(wwtp_data.ndn_plants)
    wwtp_data.curate_by_target_variable()
    wwtp_data.set_id_name('Combined_ID')
    wwtp_data.reduce_for_modelling(no_curation=False,
                                   avoid_sorbing_and_volatile=False, only_enough_data=True,
                                   only_above_LOQ=False, no_formation=False, avoid_high_std=True)
    wwtp_data.select_modeling_data()
    wwtp_data.create_modelling_input()
    # # wwtp_data.load_data('model_data')
    # # wwtp_data.randomize_y()
    #
    # # ------------------------------------#
    # # Descriptor calculation              #
    # # ------------------------------------#
    # # calculate descriptors
    descriptors = Descriptors(pep)
    descriptors.set_data(wwtp_data)
    descriptors.load_descriptors(from_csv=False, MACCS=True, mordred=True)
    # # descriptors.generate_pseudodescriptors()
    #
    # ------------------------------------#
    # Modeling                            #
    # ------------------------------------#

    wwtp_modeling = Modeling(pep, wwtp_data, descriptors)
    # wwtp_modeling.nested_cross_val_screening(feature_space_list=['maccs',
    #                                                               'maccs+ep_trig',
    #                                                               'mordred',
    #                                                               'all'],
    #                                                                config='default_wwtp')
    desired_wwtp_model = wwtp_modeling.build_final_model(regressor_name='RF',
                                                         feature_space='maccs', config='wwtp_optimized')
    print('done')
