import pathlib
from os import path

class Pepper():
    def __init__(self):
        self.root_directory = str(pathlib.Path(path.expanduser("~")))
        self.data_directory = self.root_directory + '/pepper_data/'
        self.build_directory_structure()
        self.enviPath_python = self.root_directory + '/enviPath-python/'
        self.tag = ''
        self.data_type = ''

        # todo
        self.target_variable_name = 'mean'
        self.target_variable_std_name = 'std'
        self.smiles_name = 'SMILES'
        self.id_name = 'ID'
        self.compound_name = 'compound_name'

    def set_root_directory(self, root_directory: str):
        self.root_directory = root_directory
        self.data_directory = self.root_directory + self.data_directory
        self.build_directory_structure()

    def set_data_directory(self, data_directory: str):
        self.data_directory = self.root_directory + data_directory
        self.build_directory_structure()

    def build_directory_structure(self):
        pathlib.Path(self.data_directory).mkdir(parents=True, exist_ok=True)

    def get_data_directory(self):
        return self.data_directory

    def get_root_directory(self):
        return self.root_directory

    def set_path_to_enviPath_python(self, path):
        self.enviPath_python = path

    def get_path_to_enviPath_python(self):
        return self.enviPath_python

    def set_data_type(self, data_type):
        self.data_type = data_type

    def get_data_type(self):
        return self.data_type

    def set_tag(self, tag):
        self.tag = tag

    def get_tag(self):
        return self.tag

    def set_target_variable_name(self, target_variable):
        self.set_target_variable_name = target_variable

    def set_target_variable_std_name(self, target_std_variable):
        self.set_target_variable_std_name = target_std_variable

    def set_smiles_name(self, smiles):
        self.smiles_name = smiles

    def set_id_name(self, id):
        self.id_name = id

    def build_output_filename(self, filename_string, suffix = '.tsv'):
        """
        Create filename based on current directory, data typ, and tag
        :param filename_string: file name as a string
        :return:
        """
        file_string = r'{}' + filename_string + '_{}_{}' + suffix
        complete_filename = file_string.format(self.data_directory, self.data_type, self.tag)
        return complete_filename

