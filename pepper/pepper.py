import pathlib
import os
import re

class Pepper:
    def __init__(self, renku=False):
        if renku:
            self.root_directory = str(pathlib.Path(__file__).parent.resolve()) + '/../'
            self.data_directory = os.path.join(self.root_directory, 'output')
        else:
            self.root_directory = str(pathlib.Path(os.path.expanduser("~")))
            self.data_directory = os.path.join(self.root_directory, 'pepper_data')
        self.build_directory_structure()
        self.tag = 'my_data_tag' # user-defined tag, e.g., test_data, all_data, curated_data
        self.data_type = 'unspecified_data_type' # soil, sediment, sludge, WWTP etc.
        self.setup_name = 'default_setup' # string to distinguish between different settings used, versatile
        self.curation_type = 'default_curation'

        self.target_variable_name = 'mean'
        self.target_variable_std_name = ''
        self.smiles_name = 'SMILES'
        self.id_name = 'ID'
        self.compound_name = 'Name'

        self.random_state = 42

        # development only
        # normally, enviPath-python API is installed via pip, but local git repo can be used for development purpose
        self.enviPath_python_git = os.path.join(self.root_directory, 'enviPath-python')

    def get_object_name(self):
        return re.findall(r"\.([A-Za-z]*)'", str(type(self)))[0]

    def set_root_directory(self, root_directory: str):
        self.root_directory = root_directory
        self.build_directory_structure()

    def set_data_directory(self, data_directory: str):
        self.data_directory = os.path.join(self.root_directory, data_directory)
        self.build_directory_structure()

    def build_directory_structure(self, directory: str = ''):
        pathlib.Path(os.path.join(self.data_directory, directory)).mkdir(parents=True, exist_ok=True)

    def get_data_directory(self):
        return self.data_directory

    def get_root_directory(self):
        return self.root_directory

    def set_path_to_enviPath_python(self, my_path):
        self.enviPath_python_git = my_path

    def get_path_to_enviPath_python(self):
        return self.enviPath_python_git

    def set_data_type(self, data_type):
        """
        Set a name to the context of the data used (i.e., soil, sediment, wastewater treatment plant)
        :param data_type:
        """
        self.data_type = data_type

    def set_setup_name(self, name):
        self.setup_name = name

    def get_setup_name(self):
        return self.setup_name

    def set_curation_type(self, curation_type):
        self.curation_type = curation_type

    def get_curation_type(self):
        return self.curation_type

    def get_data_type(self):
        return self.data_type

    def set_tag(self, tag):
        """
        Set a tag to identify the dataset in use.
        :param tag:
        """
        self.tag = tag

    def get_tag(self):
        return self.tag

    def get_random_state(self):
        return self.random_state

    def set_random_state(self, state):
        assert type(state) == int, 'Random state needs to be an integer'
        self.random_state = state

    def set_target_variable_name(self, target_variable):
        """
        Set the name of the variable that one desires to predict.
        :param target_variable:
        """
        self.target_variable_name = target_variable

    def get_target_variable_name(self):
        return self.target_variable_name

    def set_target_variable_std_name(self, target_std_variable):
        self.target_variable_std_name = target_std_variable

    def has_target_variable_std(self):
        if self.target_variable_std_name == '':
            return False
        else:
            return True

    def get_target_variable_std_name(self):
        return self.target_variable_std_name

    def set_smiles_name(self, smiles):
        """
        Name of the SMILES string column
        :param smiles:
        """
        self.smiles_name = smiles

    def get_smiles_name(self):
        return self.smiles_name

    def set_compound_name(self, name):
        """
        Name of the column that contains the names of the target compounds.
        :param name:
        """
        self.compound_name = name

    def get_compound_name(self):
        return self.compound_name

    def set_id_name(self, my_id):
        self.id_name = my_id

    def get_id_name(self):
        return self.id_name

    def build_output_filename(self, filename_string, suffix='.tsv'):
        """
        Create filename based on current directory, data typ, and tag
        :param filename_string: file name as a string
        :param suffix: file type, default is tab separated values
        :return:
        """
        file_string = (filename_string + '_{}_{}' + suffix).format(self.data_type, self.tag)
        complete_filename = os.path.join(self.data_directory, file_string)
        return complete_filename
