import os.path
import pandas as pd
from pepper.pepper import Pepper
from pepper.metadata import *
from pepper.util import *
from enviPath_python import enviPath
from enviPath_python.objects import Package, Scenario

class DataStructure(Pepper):
    def __init__(self, pep: Pepper):
        super().__init__()
        self.set_data_directory('/pepper_data/data_structure/')
        self.raw_data = pd.DataFrame()      # original data as loaded from source
        self.full_data = pd.DataFrame()     # current, annotated data
        self.cpd_data = pd.DataFrame()      # reduced data with extended information
        self.model_data = pd.DataFrame()    # data used for modelling with only with standardized column names

        self.tag = pep.get_tag()
        self.pepper = pep

        # enviPATH
        self.instance_host = "https://envipath.org"
        self.envipath_package = ''

    def set_target_variable_name(self, name_of_target_variable):
        self.target_variable = name_of_target_variable #name of the column

    def get_target_variable(self):
        return self.target_variable  #name of the column

    def set_smiles_name(self, smiles_column):
        self.smiles_name = smiles_column #name of the column

    def set_index(self, index_column):
        self.index_name = index_column #name of the column

    def set_envipath_package(self, uri):
        self.envipath_package = uri

    def get_target_variable(self):
        return self.full_data[self.target_variable_name]

    def get_smiles(self):
        return self.full_data[self.smiles_name]

    def get_index(self):
        return self.full_data[self.index_name]

    def load_raw_data(self, from_csv = False):
        if from_csv:
            assert os.path.exists(self.raw_data_tsv), "Error: file {} does not exist".format(self.raw_data_tsv)
            self.raw_data = pd.read_csv(self.raw_data_tsv, sep='\t')
        else:
            self.load_raw_data_from_enviPath()

    def get_model_data(self):
        return self.model_data

    @staticmethod
    def check_for_kinetics(addinfo):
        try:
            addinfo.get_halflife().get_value()
        except AttributeError:
            return False
        else:
            return True

    # functions
    def load_raw_data_from_enviPath(self):
        eP = enviPath(self.instance_host)
        pkg = Package(eP.requester, id = self.envipath_package)
        pathways = pkg.get_pathways()
        for path in pathways:
            print(pathways.index(path), path.get_id())
            for node in path.get_nodes():
                try:
                    scenarios = node.get_scenarios()
                except:
                    pass
                for scenario in scenarios:
                    print(scenario.get_id())
                    full_scenario = Scenario(eP.requester, id=scenario.get_id())
                    add_info = full_scenario.get_additional_information()
                    if self.check_for_kinetics(add_info):
                        # load all necessary data form enviPath
                        compound = CompoundStructure(eP.requester, id=node.get_default_structure().get_id())
                        metadata = Metadata(add_info)
                        try:
                            spike_smiles = CompoundStructure(eP.requester, id=add_info.get_spikecompound().get_compoundLink()).get_smiles()
                        except:
                            spike_smiles = ''

                        self.data_dict = metadata.get_scenario_information(self.data_dict, scenario,
                                                                           compound, self.data_type, spike_smiles)

        # set target variables
        self.set_target_variable_name('reported_DT50')
        self.set_smiles_name('smiles')
        # save data
        self.raw_data = pd.DataFrame(self.data_dict)
        self.raw_data.to_csv(self.raw_data_tsv, sep='\t', index=False)

    def curate_smiles(self, from_csv = False):
        print("\n############# SMILES curation ############# ")
        if from_csv:
            self.full_data = pd.read_csv(self.full_data_tsv, sep='\t')
            print("Existing file loaded from {}".format(self.full_data_tsv))
            return
        cropped_smiles, is_composite = self.get_cropped_smiles()
        self.full_data['cropped_SMILES'] = cropped_smiles
        self.full_data['is_composite'] = is_composite
        self.full_data['canonical_SMILES'] = self.get_canonicalized_smiles()
        self.full_data['cropped_canonical_SMILES'] = self.get_cropped_canonicalize_smiles()
        self.full_data['cropped_canonical_SMILES_no_stereo']= self.get_cropped_canonicalize_smiles_no_stereo()

        self.full_data.rename(columns={'smiles': 'original_SMILES'})
        self.full_data[self.smiles_name] = self.full_data['cropped_canonical_SMILES_no_stereo'] #set smiles to work with

        self.full_data.to_csv(self.full_data_tsv, sep='\t', index=False)


# todo --> make own file
class DataStructureWWTP(DataStructure):
    def load_data(self):
        # todo
        self.set_target_variable('logk')
        self.set_smiles()
        self.save_to_csv()

