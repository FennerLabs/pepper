import numpy as np
import re
import html
from pepper_lab.util import Util

class Metadata:
    """
    This class handles additional information from enviPath scenarios
    """
    def __init__(self, additional_info, description):
        self.info = additional_info
        self.des = description  # to retrieve information on high or low organic carbon scenario

    @staticmethod  # need to be fixed
    def range_to_average(input_string):
        if type(input_string) == float:  # in case we get NaN here
            return input_string
        elif input_string == ' - ' or input_string == 'NA' or input_string == '' or input_string == ' ':
            return np.NaN
        elif ';' in input_string:
            min = float(input_string.split(';')[0])
            max = float(input_string.split(';')[1])
        else:
            min = float(input_string.split(' - ')[0])
            max = float(input_string.split(' - ')[1])
        avg = np.average([min, max])
        return avg

    @staticmethod
    def is_censored(input_string):
        if '<' in input_string:
            clean_string = input_string.replace('<', '')  # replaces all '<' with ''
            return clean_string, '<'
        if '< ' in input_string:
            clean_string = input_string.replace('< ', '')  # replaces all '< ' with ''
            return clean_string, '< '
        if '>' in input_string:
            clean_string = input_string.replace('>', '')  # replaces all '>' with ''
            return clean_string, '>'
        if '> ' in input_string:
            clean_string = input_string.replace('> ', '')  # replaces all '> ' with ''
            return clean_string, '> '
        if ',' in input_string:
            clean_string = input_string.replace(',', '.')  # replaces all ',' with '.'
            return clean_string, ','
        if '(' in input_string:
            clean_string = input_string.replace('(', '')  # replaces all '(' with ''
            return clean_string, ''
        if ')' in input_string:
            clean_string = input_string.replace(')', '')  # replaces all ')' with ''
            return clean_string, ''
        else:
            return input_string, ''

    @staticmethod
    def initiate_soil_dictionary():
        D = {'compound_id': [], 'compound_name': [], 'pathway_name': [], 'node_depth': [], 'smiles': [], 'minor_major': [], # compound
             'scenario_id': [], 'study_name': [], 'study_description': [], # study/scenario
             'halflife_raw': [], 'halflife_unit': [], 'halflife_model': [], 'halflife_comment': [], # DT50
             'spike_compound': [],
             # additional information
             'acidity': [], 'acidity_unit': [],
             'temperature': [], 'temperature_unit': [],
             'CEC': [], 'CEC_unit': [],
             'OC': [],
             'biomass_start': [], 'biomass_end': [], 'biomass': [],
             'wst_value': [],
             'wst_type': [],
             'humidity': [], 'humidity_conditions': [],
             'soil_texture': [], 'sand': [], 'silt': [], 'clay': []}
        return D

    @staticmethod
    def initiate_sludge_dictionary():
        D = {
            "compound_id": [], "compound_name": [], "pathway_name": [], 'node_depth': [], "smiles": [], # compound
            "scenario_id": [], 'study_name': [], 'study_description': [], # study/scenario
            "halflife_raw": [], "halflife_unit": [], "halflife_model_TF": [], "halflife_comment": [], "halflife_model": [], # DT50
            "rateconstant": [], "rateconstant_unit": [], "rateconstant_comment": [],  # k
            # additional information
            "acidity": [], "acidity_unit": [],
            "temperature": [], "temperature_unit": [],
            "original_sludge_amount": [], "original_sludge_amount_unit": [],
            "sludge_retention_time": [], "sludge_retention_time_unit": [], "sludge_retention_time_type": [],
            "total_suspended_solids_concentration_start": [], "total_suspended_solids_concentration_end": [],
            "total_suspended_solids_concentration_unit": [],
            "addition_of_nutrients": [], "biological_treatment_technology": [],
            "bioreactor_type": [], "bioreactor_value": [], "bioreactor_value_unit": [],
            "nitrogen_content_type": [], "nitrogen_content_influent": [],
            "oxygen_demand_type": [], "oxygen_demand_value": [],
            "oxygen_uptake_rate": [], "oxygen_uptake_rate_unit": [],
            "phosphorus_content": [],
            "redox": [],
            "source_of_liquid_matrix": [],
            "type_of_addition": [],
            "type_of_aeration": [],
            "inoculum_source": [],
            "location": [],
            "purpose_of_wwtp": [],
        }
        return D

    @staticmethod # mention the relevant data-types
    def initiate_sediment_dictionary():
            D = {'compound_id': [], 'compound_name': [], 'pathway_name': [], 'node_depth': [], 'smiles': [], 'scenario_id': [],
                'study_name': [], 'major_minor': [],
                'DT50_water': [], 'DT50_sediment': [], 'DT50_total_system': [],
                'DT50_water_comment': [], 'DT50_sediment_comment': [], 'DT50_total_system_comment': [],
                'halflife_model': [], 'halflife_comment': [], 'halflife_fit': [], 'halflife_source':[],
                'spike_compound': [],
                'acidity_water': [], 'acidity_sediment': [], 'acidity_method': [],
                'bulk_density': [], 'bulk_density_unit': [],
                'column_height_water': [], 'column_height_sediment': [],
                'oxygen_content_water_start': [], 'oxygen_content_water_end': [],
                'oxygen_content_water': [],  # avg oxygen water content in water
                'oxygen_content_water_unit': [],
                'oxygen_content_sediment_start': [], 'oxygen_content_sediment_end': [],
                'oxygen_content_sediment': [],  # no oxygen content sediment value in current dataset
                'oxygen_content_sediment_unit': [],
                'CEC': [],
                'study_description': [],
                'OC_1': [], 'OC_2': [],  # range of Organic Carbon (OC) in sediment
                'OC': [],  # OC avg (average of OC_1 and OC_2), as some erroneous values on website
                'OC_type': [],
                'OM_1': [], 'OM_2': [],  # range of Organic Matter (OM) in sediment
                'OM': [],  # OM avg (average of OM_1 and OM_2)
                'TOC_1': [], 'TOC_2': [],  # range of Total Organic carbon (TOC) in water layer
                'TOC': [],  # TOC avg (average of TOC_1 and TOC_2)
                'TOC_unit': [],
                'DOC_1': [], 'DOC_2': [],  # range of Dissolved Organic carbon (DOC) in water layer
                'DOC': [],  # DOC avg (average of DOC_1 and DOC_2)
                'DOC_unit': [],
                'redox_water_start': [], 'redox_water_end': [],  # In water, redox potential start and end, respectively
                # In sediment, redox potential start and end, respectively. Only few values in the current dataset
                'redox_sediment_start': [], 'redox_sediment_end': [],
                'biomass_cells_count_water': [],
                'biomass_cells_count_water_unit': [],
                'biomass_cells_count_sediment': [],
                'biomass_cells_count_sediment_unit': [],
                'biomass_sediment_mg_start': [], 'biomass_sediment_mg_end': [],  # range of biomass at start and end in sediment
                'biomass': [],  # avg biomass (average of biomass_sediment_start and biomass_sediment_end)
                'biomass_sediment_unit': [],
                'temperature': [], 'temperature_unit': [],
                'sample_location': [],
                'sediment_porosity': [],
                'initial_sediment_mass': [],
                'initial_sediment_mass_wet_or_dry': [],
                'sediment_condition': [],
                'initial_volume_water': [], 'initial_volume_water_unit': [],
                'soil_texture': [], 'sand': [], 'silt': [], 'clay': []}
            return D



    def fetch_mean_value(self, name, key1, key2):
        try:
            raw = self.info[name].params
            val1 = float(raw[key1])
            val2 = float(raw[key2])
        except:
            return np.NaN
        else:
            if np.isnan(val1) and np.isnan(val2):
                return np.nan
            elif np.isnan(val1):
                mean = val2
            elif np.isnan(val2):
                mean = val1
            else:
                mean = np.average([val1, val2])
            return np.round(mean, 4)

    def fetch_normal_value(self, name, key, typ):
        default_return = "" if typ is str else np.nan
        try:
            value = typ(self.info[name].params[key])
        except:
            return default_return
        else:
            if key == "unit" or "type" in key:
                return html.unescape(value)
            return value

    # Fetch bulk density in kg/m3
    def fetch_bulk_density(self):
        try:
            bulk_density = self.info.get_bulkdens().get_value()
        except:
            return np.NaN
        else:
            return bulk_density

    def fetch_bulk_density_unit(self):
        try:
            unit = self.info.get_bulkdens().get_unit()
        except:
            return ''
        else:
            return unit


    def fetch_cec(self):
        try:
            cec = self.info.get_cec().get_value()
        except:
            return np.NaN
        else:
            # if ',' in cec:  # 2-3 values had a comma in place of decimal point
            # but this cannot be used as cec data type is float
            # hence, those values were corrected manually in the dataset
            #     cec = cec.replace(',', '.')
            if cec == '2023-06-14 00:00:00':  # todo: remove once website database is fixed
                cec = 1.46  # This value taken from nearest scenario of the compound (Foramsulfuron)
                # and it needs to be verified from its DAR
            # cec, var = self.is_censored(cec)
            return cec

    def fetch_organic_content(self):
        try:
            raw = self.info.get_omcontent().get_value()
        except:
            return np.NaN
        else:
            raw_list = raw.split(';')
            oc = np.NaN
            for i in raw_list:
                if i == 'OC':
                    oc = val
                elif i == 'OM':
                    oc = val / 1.7  # OC = OM / 1.7, source: Schwarzenbach
                else:
                    if '<' in i:
                        val = float(i[1:])
                        print("Warning: {} was converted to {}".format(i, val))
                    elif i == '' or i == '-':
                        val = np.NaN
                    else:
                        val = float(i)
            return oc


    def fetch_biomass(self):
        try:
            raw = self.info.get_biomass().get_value()
        except:
            return np.NaN, np.NaN
        else:
            l = raw.split(' - ')
            return float(l[0]), float(l[1])


    def fetch_temperature(self):
        try:
            raw = self.info.get_temperature().get_value()
        except:
            return np.NaN
        else:
            if raw == ' ' or raw == '' or raw == ';' or raw == ' ;' or raw == '; ':
                return np.NaN
            else:
                # min = float(raw.split(';')[0])
                # max = float(raw.split(';')[1])
                min = raw.split(';')[0]
                max = raw.split(';')[1]
                if min == ' ' or min == '':
                    min = np.NaN
                    return min
                if '-' in min:
                    min = min.split('-')[0]
                    return min
                if '-' in max:
                    max = max.split('-')[1]
                    return max
                if max == ' ' or max == '':
                    max = np.NaN
                    return max
                if ',' in min or ',' in max:
                    min = min.replace(',', '.')
                    max = max.replace(',', '.')
            return np.round(np.average([float(min), float(max)]), 0)

    def fetch_wst(self):
        try:
            raw = self.info.get_waterstoragecapacity().get_value()
        except:
            return np.NaN, ''
        else:
            raw_list = raw.replace(" ", "").split('-')
            if len(raw_list) < 4:
                value = float(raw_list[0])
                type = raw_list[1]
            else:
                value = np.NaN
                type = raw_list[2]
            return value, type

    def fetch_humidity(self):
        try:
            raw = self.info.get_humidity().get_value()
        except:
            return np.NaN, ''
        else:
            if type(raw) == float:
                return raw, ''
            else:
                l = raw.split(' - ')
                return float(l[0]), l[1]


    # fetch sample location
    def fetch_sample_location(self):
        try:
            location = self.info.get_samplelocation().get_value()
            # map plot function
        except:
            return ''
        else:
            return location

    def fetch_soiltexture1(self):
        try:
            raw = self.info.get_soiltexture1().get_value()
        except:
            return ''
        else:
            return raw

    def fetch_soiltexture2(self):
        try:
            raw = self.info.get_soiltexture2().get_value()
        except:
            return np.NaN, np.NaN, np.NaN
        else:
            values = re.findall(r'\s([\d.]+)%', raw)  #or values =
            if values == []:
                return np.NaN, np.NaN, np.NaN
            elif len(values) < 3 and 'E' in raw:
                values = re.findall(r'\s([E\-\d.]+)%', raw)
                return self.get_float_or_nan(values[0]), self.get_float_or_nan(values[1]), self.get_float_or_nan(values[2])
            elif '<' in raw:
                new_val = self.is_censored(raw)
                values = new_val.split(';')
                return self.get_float_or_nan(values[0]), self.get_float_or_nan(values[1]), self.get_float_or_nan(
                    values[2])
            else:
                return self.get_float_or_nan(values[0]), self.get_float_or_nan(values[1]), self.get_float_or_nan(
                    values[2])
              # sand, silt, clay





    @staticmethod
    def get_float_or_nan(x):
        try:
            return float(x)
        except:
            return np.NaN

    def get_scenario_information(self, D, scenario, compound, data_type, spike_smiles, description, path_name, depth):
        if data_type == 'soil':
            D = self.get_soil_scenario_information(D, scenario, compound, spike_smiles, description, path_name, depth)
        elif data_type == 'sediment':
            D = self.get_sediment_scenario_information(D, scenario, compound, spike_smiles, description, path_name, depth)
        elif data_type == 'sludge':
            D = self.get_sludge_scenario_information(D, scenario, compound, description, path_name, depth)
        else:
            raise NotImplementedError
        return D

    def get_sludge_scenario_information(self, D, scenario, compound, description, path_name, depth):
        if D == {}:
            D = self.initiate_sludge_dictionary()
        # Compound informatin
        D['compound_id'].append(compound.get_id())
        D['compound_name'].append(compound.get_name())

        D['smiles'].append(compound.get_smiles())
        # Scenario/study information
        D['scenario_id'].append(scenario.get_id())
        D['study_name'].append(scenario.get_name().split(' - ')[0])
        D['pathway_name'].append(path_name)
        D['node_depth'].append(depth)

        D['study_description'].append(description)
        if "acidity" in self.info:
            D['acidity'].append(self.fetch_mean_value("acidity", 'lowPh', 'highPh'))
            D['acidity_unit'].append(self.info['acidity'].get_unit())
        else:
            D['acidity'].append(np.NaN)
            D['acidity_unit'].append(np.NaN)
        if "additionofnutrients" in self.info:
            D['addition_of_nutrients'].append(self.info['additionofnutrients'].get_additionofnutrients())
        else:
            D['addition_of_nutrients'].append(np.nan)
        if "biologicaltreatmenttechnology" in self.info:
            D['biological_treatment_technology'].append(self.info['biologicaltreatmenttechnology'].get_biologicaltreatmenttechnology())
        else:
            D['biological_treatment_technology'].append(np.nan)

        if "bioreactor" in self.info:
            D['bioreactor_type'].append(self.info['bioreactor'].get_bioreactortype())
            D['bioreactor_value'].append(self.info['bioreactor'].get_bioreactorsize())
            D['bioreactor_value_unit'].append(self.info['bioreactor'].get_unit())
        else:
            D['bioreactor_type'].append(np.nan)
            D['bioreactor_value'].append(np.nan)
            D['bioreactor_value_unit'].append(np.nan)
        if "halflife" in self.info:
            D['halflife_raw'].append(np.nanmean([self.info["halflife"].get_lower(), self.info["halflife"].get_upper()]))
            D['halflife_unit'].append(self.info['halflife'].get_unit())
            D['halflife_comment'].append(self.info['halflife'].get_comment())
            if self.info['halflife'].get_firstOrder():
                D['halflife_model_TF'].append('SFO')
            else:
                D['halflife_model_TF'].append('nan')
        else:
            D['halflife_raw'].append(np.nan)
            D['halflife_unit'].append(np.nan)
            D['halflife_comment'].append(np.nan)
            D['halflife_model_TF'].append(np.nan)
        if "inoculumsource" in self.info:
            D['inoculum_source'].append(self.info['inoculumsource'].get_inoculumsource())
        else:
            D['inoculum_source'].append(np.nan)
        if "location" in self.info:
            D['location'].append(self.info['location'].get_location())
        else:
            D['location'].append(np.nan)
        if "nitrogencontent" in self.info:
            D['nitrogen_content_influent'].append(self.info['nitrogencontent'].get_nitrogencontentInfluent())
            D['nitrogen_content_type'].append(self.info['nitrogencontent'].get_nitrogencontentType())
        else:
            D['nitrogen_content_influent'].append(np.nan)
            D['nitrogen_content_type'].append(np.nan)
        if "originalsludgeamount" in self.info:
            D['original_sludge_amount'].append(self.info['originalsludgeamount'].get_originalsludgeamount())
            D['original_sludge_amount_unit'].append(self.info['originalsludgeamount'].get_unit())
        else:
            D['original_sludge_amount'].append(np.nan)
            D['original_sludge_amount_unit'].append(np.nan)
        if "oxygendemand" in self.info:
            D['oxygen_demand_type'].append(self.info['oxygendemand'].get_oxygendemandType())
            D['oxygen_demand_value'].append(self.info['oxygendemand'].get_oxygendemandInfluent())
        else:
            D['oxygen_demand_type'].append(np.nan)
            D['oxygen_demand_value'].append(np.nan)
        if "oxygenuptakerate" in self.info:
            start = self.info['oxygenuptakerate'].get_oxygenuptakerateStart()
            end = self.info['oxygenuptakerate'].get_oxygenuptakerateEnd()
            D['oxygen_uptake_rate_unit'].append(self.info['oxygenuptakerate'].get_unit())
            if start and end:

                D['oxygen_uptake_rate'].append(np.nanmean([float(start), float(end)]))
            elif start and not end:
                D['oxygen_uptake_rate'].append(float(self.info['oxygenuptakerate'].get_oxygenuptakerateStart()))
            else:
                D['oxygen_uptake_rate'].append(float(self.info['oxygenuptakerate'].get_oxygenuptakerateEnd()))
        else:
            D['oxygen_uptake_rate_unit'].append(np.nan)
            D['oxygen_uptake_rate'].append(np.nan)
        if "phosphoruscontent" in self.info:
            D['phosphorus_content'].append(self.info['phosphoruscontent'].get_phosphoruscontentInfluent())
        else:
            D['phosphorus_content'].append(np.nan)
        if "purposeofwwtp" in self.info:
            D['purpose_of_wwtp'].append(self.info['purposeofwwtp'].get_purposeofwwtp())
        else:
            D['purpose_of_wwtp'].append(np.nan)
        if "rateconstant" in self.info:
            D['rateconstant'].append(np.nanmean([self.info['rateconstant'].get_rateconstantlower(), self.info['rateconstant'].get_rateconstantupper()]))
            D['rateconstant_unit'].append(self.info['rateconstant'].get_unit())
            D['halflife_model'].append(self.info['rateconstant'].get_rateconstantorder())
            D['rateconstant_comment'].append(self.info['rateconstant'].get_rateconstantcomment())
        else:
            D['rateconstant'].append(np.nan)
            D['rateconstant_unit'].append(np.nan)
            D['halflife_model'].append(np.nan)
            D['rateconstant_comment'].append(np.nan)
        if "redox" in self.info:
            D['redox'].append(self.info['redox'].get_redoxType())
        else:
            D['redox'].append(np.nan)
        if "sludgeretentiontime" in self.info:
            D['sludge_retention_time'].append(self.info['sludgeretentiontime'].get_sludgeretentiontime())
            D["sludge_retention_time_unit"].append(self.info['sludgeretentiontime'].get_unit())
            D['sludge_retention_time_type'].append(self.info['sludgeretentiontime'].get_sludgeretentiontimeType())
        else:
            D['sludge_retention_time'].append(np.nan)
            D['sludge_retention_time_unit'].append(np.nan)
            D['sludge_retention_time_type'].append(np.nan)
        if "sourceofliquidmatrix" in self.info:
            D['source_of_liquid_matrix'].append(self.info['sourceofliquidmatrix'].get_sourceofliquidmatrix())
        else:
            D['source_of_liquid_matrix'].append(np.nan)
        if "temperature" in self.info:
            D['temperature'].append(np.nanmean([float(self.info['temperature'].get_temperatureMin()), float(self.info['temperature'].get_temperatureMax())]))
            D['temperature_unit'].append(self.info['temperature'].get_unit())
        else:
            D['temperature'].append(np.nan)
            D['temperature_unit'].append(np.nan)
        if "tts" in self.info:
            D['total_suspended_solids_concentration_start'].append(self.info['tts'].get_ttsStart())
            D['total_suspended_solids_concentration_end'].append(self.info['tts'].get_ttsEnd())
            D['total_suspended_solids_concentration_unit'].append(self.info['tts'].get_unit())
        else:
            D['total_suspended_solids_concentration_start'].append(np.nan)
            D['total_suspended_solids_concentration_end'].append(np.nan)
            D['total_suspended_solids_concentration_unit'].append(np.nan)
        if "typeofaddition" in self.info:
            D['type_of_addition'].append(self.info['typeofaddition'].get_typeofaddition())
        else:
            D['type_of_addition'].append(np.nan)
        if "typeofaeration" in self.info:
            D['type_of_aeration'].append(self.info['typeofaeration'].get_typeofaeration())
        else:
            D['type_of_aeration'].append(np.nan)
        return D

    def get_soil_scenario_information(self, D, scenario, compound, spike_smiles, description, path_name,depth):
        if D == {}:
            D = self.initiate_soil_dictionary()
        # compound info
        D['spike_compound'].append(spike_smiles)
        D['compound_id'].append(compound.get_id())
        D['compound_name'].append(compound.get_name())
        D['smiles'].append(compound.get_smiles())
        if 'minormajor' in self.info:
            D['minor_major'].append(self.info['minormajor'].get_radiomin())
        else:
            D['minor_major'].append(np.nan)
        # D['spike_compound'].append(spike_smiles)
        # study
        D['scenario_id'].append(scenario.get_id())
        D['study_name'].append(scenario.get_name().split(' - ')[0])
        D['study_description'].append(description)
        D['pathway_name'].append(path_name)
        D['node_depth'].append(depth)
        # add halflife details
        D['halflife_raw'].append(self.fetch_mean_value("halflife", "lower", "upper")) # self.info['halflife']
        if 'halflife' in self.info:
            D['halflife_unit'].append(self.info['halflife'].get_unit())
            if self.info['halflife'].get_firstOrder():
                D['halflife_model'].append('SFO')
            else:
                D['halflife_model'].append('nan')
            D['halflife_comment'].append(self.info['halflife'].get_comment())
        
        else:
            D['halflife_unit'].append(np.nan)
            D['halflife_model'].append(np.nan)
            D['halflife_comment'].append(np.nan)
        # fetch additional information
        if "acidity" in self.info:
            D['acidity'].append(self.fetch_mean_value("acidity", 'lowPh', 'highPh'))
            D['acidity_unit'].append(self.info['acidity'].get_unit())
        else:
            D['acidity'].append(np.NaN)
            D['acidity_unit'].append(np.NaN)

        if "temperature" in self.info:
            D['temperature'].append(self.fetch_mean_value("temperature", "temperatureMin", "temperatureMax"))
            D['temperature_unit'].append(self.info['temperature'].get_unit())
        else:
            D['temperature'].append(np.NaN)
            D['temperature_unit'].append(np.NaN)

        if "cec" in self.info:
            D['CEC'].append(self.info['cec'].get_cecdata())  # cation exchange capacity
            D['CEC_unit'].append(self.info['cec'].get_unit())

        else:
            D['CEC'].append(np.NaN)
            D['CEC_unit'].append(np.NaN)
        if "organiccontent" in self.info:
            D['OC'].append(self.fetch_mean_value('organiccontent','OC_content_low', 'OC_content_high'))  
        
        else:
            D['OC'].append(np.NaN)
        if "biomass" in self.info:
            D['biomass_start'].append(self.info['biomass'].get_biomassStart())
            D['biomass_end'].append(self.info['biomass'].get_biomassEnd())
            D['biomass'].append(np.round(np.average([D['biomass_start'], D['biomass_end']]), 2))
        else:
            D['biomass_start'].append(np.NaN)
            D['biomass_end'].append(np.NaN)
            D['biomass'].append(np.NaN)
        if "waterstoragecapacity" in self.info:
            D['wst_value'].append(self.info['waterstoragecapacity'].get_maximumWaterstoragecapacity())
            D['wst_type'].append(self.info['waterstoragecapacity'].get_wstConditions())
        else:
            D['wst_value'].append(np.NaN)
            D['wst_type'].append(np.NaN)
        if "humidity" in self.info:
            D['humidity'].append(self.info['humidity'].get_expHumid())
            D['humidity_conditions'].append(self.info['humidity'].get_humConditions())
        else:
            D['humidity'].append(np.NaN)
            D['humidity_conditions'].append(np.NaN)
        if "soiltexture1" in self.info:
            D['soil_texture'].append(self.info['soiltexture1'].get_soilTextureType())
        else:
            D['soil_texture'].append(np.NaN)
        if "soiltexture2" in self.info:
            D['sand'].append(self.info['soiltexture2'].get_sand())
            D['silt'].append(self.info['soiltexture2'].get_silt())
            D['clay'].append(self.info['soiltexture2'].get_clay())
        else:
            D['sand'].append(np.NaN)
            D['silt'].append(np.NaN)
            D['clay'].append(np.NaN)
        return D

    def get_sediment_scenario_information(self, D, scenario, compound, spike_smiles, description, path_name,depth):
        # compound info
        if D == {}:
            D = self.initiate_sediment_dictionary()
        D['compound_id'].append(compound.get_id())
        D['compound_name'].append(compound.get_name())
        D['smiles'].append(compound.get_smiles())
        D['spike_compound'].append(spike_smiles)
        D['scenario_id'].append(scenario.get_id())
        D['pathway_name'].append(path_name)
        D['node_depth'].append(depth)
        # System information: High OC or Low OC system
        D['study_description'].append(description)
        D['OC_type'].append("see description")
        # Fetch other data points

        D['study_name'].append(scenario.get_name().split(' - ')[0])
        if 'minormajor' in self.info:
            D['major_minor'].append(self.info['minormajor'].get_radiomin())
        else:
            D['major_minor'].append(np.nan)
        # fetch halflife details - total system, water, sediment, model_type, comment, fit
        if "halflife_ws" in self.info:
            D['DT50_total_system'].append(np.nanmean([self.info['halflife_ws'].get_total_high(), self.info['halflife_ws'].get_total_low()]))
            
            water_low = self.info['halflife_ws'].get_water_low()
            water_high = self.info['halflife_ws'].get_water_high()
            sediment_low = self.info['halflife_ws'].get_sediment_low()
            sediment_high = self.info['halflife_ws'].get_sediment_high()
            if water_low and water_high:
                D['DT50_water'].append(np.nanmean([water_low, water_high]))
            else:
                D['DT50_water'].append(np.nan)
            if sediment_low and sediment_high:
                D['DT50_sediment'].append(np.nanmean([sediment_low, sediment_high]))
            else:
                D['DT50_sediment'].append(np.nan)
            
            D['DT50_water_comment'].append(self.info['halflife_ws'].get_comment_ws())
            D['DT50_total_system_comment'].append(self.info['halflife_ws'].get_comment_ws())
            D['DT50_sediment_comment'].append(self.info['halflife_ws'].get_comment_ws())
            D['halflife_model'].append(self.info['halflife_ws'].get_model_ws())
            D['halflife_comment'].append(self.info['halflife_ws'].get_comment_ws())
            D['halflife_fit'].append(self.info['halflife_ws'].get_fit_ws())
            D['halflife_source'].append(self.info['halflife_ws'].get_source_ws())
        else:
            D['DT50_total_system'].append(np.nan)
            D['DT50_water'].append(np.nan)
            D['DT50_water_comment'].append(np.nan)
            D['DT50_total_system_comment'].append(np.nan)
            D['DT50_sediment_comment'].append(np.nan)
            D['DT50_sediment'].append(np.nan)
            D['halflife_model'].append(np.nan)
            D['halflife_comment'].append(np.nan)
            D['halflife_fit'].append(np.nan)
            D['halflife_source'].append(np.nan)
        #  fetch pH values for surface water and sediment, and the method used for measuring pH in sediment
        if "temperature" in self.info:
            D['temperature'].append(np.nanmean([float(self.info['temperature'].get_temperatureMin()), float(self.info['temperature'].get_temperatureMax())]))
            D['temperature_unit'].append(self.info['temperature'].get_unit())
        else:
            D['temperature'].append(np.NaN)
            D['temperature_unit'].append(np.NaN)
        if "acidity_ws" in self.info:
            if self.info['acidity_ws'].get_pH_water_low() or self.info['acidity_ws'].get_pH_water_high():
                D['acidity_water'].append(np.nanmean([self.info['acidity_ws'].get_pH_water_low(), self.info['acidity_ws'].get_pH_water_high()])) # pH surface water low carbon
            else:
                D['acidity_water'].append(np.nan)
            if self.info['acidity_ws'].get_pH_sediment_low() or self.info['acidity_ws'].get_pH_sediment_high():
                D['acidity_sediment'].append(np.nanmean([self.info['acidity_ws'].get_pH_sediment_low(),self.info['acidity_ws'].get_pH_sediment_high()]))  # pH sediment
            else:
                D['acidity_sediment'].append(np.nan)
            D['acidity_method'].append(self.info['acidity_ws'].get_acidityType())  # pH method in sediment
        else:
            D['acidity_water'].append(np.nan)
            D['acidity_sediment'].append(np.nan)
            D['acidity_method'].append(np.nan)
        if "bulkdens" in self.info:
            D['bulk_density'].append(self.info['bulkdens'].get_bulkdensity())
            D['bulk_density_unit'].append(self.info['bulkdens'].get_unit())
        else:
            D['bulk_density'].append(np.nan)
            D['bulk_density_unit'].append(np.nan)
        if 'cec' in self.info:
            D['CEC'].append(self.info['cec'].get_cecdata())
        else:
            D['CEC'].append(np.nan)
        # column height for water and sediment phases, respectively
        if "columnheight" in self.info:
            D['column_height_water'].append(self.info['columnheight'].get_column_height_water())
            D['column_height_sediment'].append(self.info['columnheight'].get_column_height_sediment())
        else:
            D['column_height_water'].append(np.nan)
            D['column_height_sediment'].append(np.nan)
        # initial sediment mass (dry/wet)
        if 'initialmasssediment' in self.info:
            D['initial_sediment_mass'].append(self.info['initialmasssediment'].get_initial_mass_sediment())
            D['initial_sediment_mass_wet_or_dry'].append(self.info['initialmasssediment'].get_wet_or_dry())
        else:
            D['initial_sediment_mass'].append(np.nan)
            D['initial_sediment_mass_wet_or_dry'].append(np.nan)
        # initial volume of water
        if 'initialvolumewater' in self.info:
            D['initial_volume_water'].append(self.info['initialvolumewater'].get_initialvolumewater())
            D['initial_volume_water_unit'].append(self.info['initialvolumewater'].get_unit())
        else:
            D['initial_volume_water'].append(np.nan)
            D['initial_volume_water_unit'].append(np.nan)
        # Organic carbon in water layer - Total Organic Carbon [TOC] values % and Dissolved Organic Carbon [DOC] values 
        if 'organiccarbonwater' in self.info:
            toc1 = self.info['organiccarbonwater'].get_TOC_low()
            toc2 = self.info['organiccarbonwater'].get_TOC_high()
            doc1 = self.info['organiccarbonwater'].get_DOC_low()
            doc2 = self.info['organiccarbonwater'].get_DOC_high()
            if toc1 or toc2:
                toc1 = toc1.replace('-', '.').replace('<', '')
                toc2 = toc2.replace('-', '.').replace('<', '')
                D['TOC_1'].append(self.info['organiccarbonwater'].get_TOC_low())
                D['TOC_2'].append(self.info['organiccarbonwater'].get_TOC_high())
                D['TOC'].append(np.nanmean([float(toc1), float(toc2)]))
                D['TOC_unit'].append(self.info['organiccarbonwater'].get_unit())
            else:
                D['TOC_1'].append(np.nan)
                D['TOC_2'].append(np.nan)
                D['TOC'].append(np.nan)
                D['TOC_unit'].append(np.nan)
            if doc1 or doc2:
                doc1 = doc1.replace('-', '.').replace('<', '')
                doc2 = doc2.replace('-', '.').replace('<', '')
                D['DOC_1'].append(self.info['organiccarbonwater'].get_DOC_low())
                D['DOC_2'].append(self.info['organiccarbonwater'].get_DOC_high())
                D['DOC'].append(np.nanmean([float(doc1), float(doc2)]))
                D['DOC_unit'].append(self.info['organiccarbonwater'].get_unit())
            else:
                D['DOC_1'].append(np.nan)
                D['DOC_2'].append(np.nan)
                D['DOC'].append(np.nan)
                D['DOC_unit'].append(np.nan)
      
        else:
            D['TOC_1'].append(np.nan)
            D['TOC_2'].append(np.nan)
            D['TOC'].append(np.nan)
            D['TOC_unit'].append(np.nan)
            D['DOC_1'].append(np.nan)
            D['DOC_2'].append(np.nan)
            D['DOC'].append(np.nan)
            D['DOC_unit'].append(np.nan)
        # Organic content in sediment organic carbon [OC] and organic matter [OM] values
        if "organiccontent" in self.info:
            oc1 = self.info["organiccontent"].get_OC_content_low()
            oc2 = self.info["organiccontent"].get_OC_content_high()
            if oc1 or oc2:
                oc1 = oc1.replace('<', '')
                oc2 = oc2.replace('<', '')
                D['OC_1'].append(self.info["organiccontent"].get_OC_content_low())
                D['OC_2'].append(self.info["organiccontent"].get_OC_content_high())
                D['OC'].append(np.nanmean([float(oc1), float(oc2)]))
            else:
                D['OC_1'].append(np.nan)
                D['OC_2'].append(np.nan)
                D['OC'].append(np.nan)
            if self.info["organiccontent"].get_OM_content_low() or self.info["organiccontent"].get_OM_content_high():
                
                D['OM_1'].append(self.info["organiccontent"].get_OM_content_low())
                D['OM_2'].append(self.info["organiccontent"].get_OM_content_high())
                D['OM'].append(np.nanmean([float(self.info["organiccontent"].get_OM_content_low()), float(self.info["organiccontent"].get_OM_content_high())]))
            else:
                D['OM_1'].append(np.nan)
                D['OM_2'].append(np.nan)
                D['OM'].append(np.nan)
        else:
            D['OC_1'].append(np.nan)
            D['OC_2'].append(np.nan)
            D['OC'].append(np.nan)
            D['OM_1'].append(np.nan)
            D['OM_2'].append(np.nan)
            D['OM'].append(np.nan)
        # Oxygen content in water layer
        if "oxygencontent" in self.info:
            ox_w_low = self.info["oxygencontent"].get_oxygen_content_water_low()
            ox_w_high = self.info["oxygencontent"].get_oxygen_content_water_high()
            ox_s_low = self.info["oxygencontent"].get_oxygen_content_sediment_low()
            ox_s_high = self.info["oxygencontent"].get_oxygen_content_sediment_high()
            if ox_w_low == 'NA':
                ox_w_low = np.nan
            if ox_w_high == 'NA':
                ox_w_high = np.nan
            if ox_s_low == 'NA':
                ox_s_low = np.nan
            if ox_s_high == 'NA':
                ox_s_high = np.nan
            if ox_w_low and ox_w_high :
        
                D["oxygen_content_water_start"].append(float(ox_w_low))
                D["oxygen_content_water_end"].append(float(ox_w_high))
                D["oxygen_content_water"].append(np.nanmean([float(ox_w_low), float(ox_w_high)]))
                D["oxygen_content_water_unit"].append(self.info["oxygencontent"].get_unit())
            else:
                D["oxygen_content_water_start"].append(np.nan)
                D["oxygen_content_water_end"].append(np.nan)
                D["oxygen_content_water"].append(np.nan)
                D["oxygen_content_water_unit"].append(np.nan)
    
            if ox_s_low and ox_s_high:      
                D["oxygen_content_sediment_start"].append(float(ox_s_low))
                D["oxygen_content_sediment_end"].append(float(ox_s_high))
                D["oxygen_content_sediment"].append(np.nanmean([float(ox_s_low), float(ox_s_high)]))
                D["oxygen_content_sediment_unit"].append(self.info["oxygencontent"].get_unit())
            else:
                D["oxygen_content_sediment_start"].append(np.nan)
                D["oxygen_content_sediment_end"].append(np.nan)
                D["oxygen_content_sediment"].append(np.nan)
                D["oxygen_content_sediment_unit"].append(np.nan)
       
        else:
            D["oxygen_content_water_start"].append(np.nan)
            D["oxygen_content_water_end"].append(np.nan)
            D["oxygen_content_water"].append(np.nan)
            D["oxygen_content_water_unit"].append(np.nan)
            D["oxygen_content_sediment_unit"].append(np.nan)
            D["oxygen_content_sediment_start"].append(np.nan)
            D["oxygen_content_sediment_end"].append(np.nan)
            D["oxygen_content_sediment"].append(np.nan)
        
        if "biomass_ws" in self.info:
            if self.info['biomass_ws'].get_start_water_cells() or self.info['biomass_ws'].get_end_water_cells():
                D['biomass_cells_count_water'].append(np.nanmean([self.info['biomass_ws'].get_start_water_cells(), self.info['biomass_ws'].get_end_water_cells()]))
                D['biomass_cells_count_water_unit'].append(self.info['biomass_ws'].get_unit())
            else: 
                D['biomass_cells_count_water'].append(np.nan)
                D['biomass_cells_count_water_unit'].append(np.nan)
            if self.info['biomass_ws'].get_start_sediment_cells() or self.info['biomass_ws'].get_end_sediment_cells():
                D['biomass_cells_count_sediment'].append(np.nanmean([float(self.info['biomass_ws'].get_start_sediment_cells()), float(self.info['biomass_ws'].get_end_sediment_cells())]))
                D['biomass_cells_count_sediment_unit'].append(self.info['biomass_ws'].get_unit())
            else:
                D['biomass_cells_count_sediment'].append(np.nan)
                D['biomass_cells_count_sediment_unit'].append(np.nan)
            start = self.info['biomass_ws'].get_start_sediment_mg()
            end = self.info['biomass_ws'].get_end_sediment_mg()
            if start and end:
                start = start.replace('<', '')
                end = end.replace('<', '')
                D['biomass_sediment_mg_start'].append(self.info['biomass_ws'].get_start_sediment_mg())
                D['biomass_sediment_mg_end'].append(self.info['biomass_ws'].get_end_sediment_mg())

                D['biomass'].append(np.nanmean([float(start), float(end)]))
                D['biomass_sediment_unit'].append("mg")
            else:
                D['biomass_sediment_mg_start'].append(np.nan)
                D['biomass_sediment_mg_end'].append(np.nan)
                D['biomass'].append(np.nan)
                D['biomass_sediment_unit'].append(np.nan)
        else:
            D['biomass_cells_count_water'].append(np.nan)
            D['biomass_cells_count_water_unit'].append(np.nan)
            D['biomass_cells_count_sediment'].append(np.nan)
            D['biomass_cells_count_sediment_unit'].append(np.nan)
            D['biomass_sediment_mg_start'].append(np.nan)
            D['biomass_sediment_mg_end'].append(np.nan)
            D['biomass'].append(np.nan)
            D['biomass_sediment_unit'].append(np.nan)
        # Redox potential of water and sediment
        if 'redoxpotential' in self.info:
            D['redox_water_start'].append(self.info['redoxpotential'].get_lowPotentialWater())
            D['redox_water_end'].append(self.info['redoxpotential'].get_highPotentialWater())
            D['redox_sediment_start'].append(self.info['redoxpotential'].get_lowPotentialSediment())
            D['redox_sediment_end'].append(self.info['redoxpotential'].get_highPotentialSediment())
        else:
            D['redox_water_start'].append(np.nan)
            D['redox_water_end'].append(np.nan)
            D['redox_sediment_start'].append(np.nan)
            D['redox_sediment_end'].append(np.nan)
        # sample location
        D['sediment_condition'].append(np.nan)
        if 'samplelocation' in self.info:
            D['sample_location'].append(self.info['samplelocation'].get_samplelocation())
        else:
            D['sample_location'].append(np.nan)
        if 'sedimentporosity' in self.info:
            D['sediment_porosity'].append(self.info['sedimentporosity'].get_sedimentporosity())
        else:
            D['sediment_porosity'].append(np.nan)
        if 'soiltexture1' in self.info:
            D['soil_texture'].append(self.info["soiltexture1"].get_soilTextureType())
        else:
            D['soil_texture'].append(np.nan)
        if 'soiltexture2' in self.info:
            D['sand'].append(self.info["soiltexture2"].get_sand())  # sand
            D['silt'].append(self.info["soiltexture2"].get_silt())  # silt                  
            D['clay'].append(self.info["soiltexture2"].get_clay())  # clay
        else:
            D['sand'].append(np.nan)
            D['silt'].append(np.nan)
            D['clay'].append(np.nan)
        return D



# these functions will return halflife in water, sediment and total system respectively; halflife model, r^2, chi^2
    def fetch_halflife_total_system_value(self):
        try:
            raw = self.info.get_halflife_ws().get_value()
        except:
            return np.NaN, ''
        else:
            if len(raw.split(';')) < 4:
                print('Warning: incomplete half-life information - {}'.format(raw))
                DT50_total_sys = np.NaN
                comment = ''
                return DT50_total_sys, comment
            else:
                DT50_total_sys = raw.split(';')[3]
                DT50_total_sys, comment = self.is_censored(DT50_total_sys)
                return DT50_total_sys, comment

    def fetch_halflife_water_value(self):
        try:
            raw = self.info.get_halflife_ws().get_value()
        except:
            return np.NaN, ''
        else:
            if len(raw.split(';')) < 4:
                print('Warning: incomplete half-life information - {}'.format(raw))
                DT50_water = np.NaN
                comment = ''
                return DT50_water, comment
            else:
                DT50_water = raw.split(';')[4]
                DT50_water, comment = self.is_censored(DT50_water)
                return DT50_water, comment

    def fetch_halflife_sediment_value(self):
        try:
            raw = self.info.get_halflife_ws().get_value()
        except:
            return np.NaN, ''
        else:
            if len(raw.split(';')) < 4:
                print('Warning: incomplete half-life information - {}'.format(raw))
                DT50_sed = np.NaN
                comment = ''
                return DT50_sed, comment
            else:
                DT50_sed = raw.split(';')[5]
                DT50_sed, comment = self.is_censored(DT50_sed)
                return DT50_sed, comment


    def fetch_halflife_ws_model(self):
        try:
            raw = self.info.get_halflife_ws().get_value()
        except:
            return ''
        else:
            return raw.split(';')[0]


    def fetch_halflife_ws_comment(self):
        try:
            raw = self.info.get_halflife_ws().get_value()

        except:
            return ''
        else:

            return raw.split(';')[2]

    def fetch_halflife_ws_fit(self):
        try:
            raw = self.info.get_halflife_ws().get_value()
        except:
            return ''
        else:
            return raw.split(';')[1]

    # fetch initial sediment mass and condition for water-sediment data
    def fetch_initial_sediment_mass(self):
        try:
            raw_sediment_mass = self.info.get_initialmasssediment().get_value()
        except:
            return np.NaN, np.NaN, ''
        else:
            initial_sediment_mass = self.get_float_or_nan(raw_sediment_mass.split(';')[0])
            sediment_condition = raw_sediment_mass.split(';')[1]
            if 'dry' in sediment_condition:
                initial_sediment_mass_dry = float(initial_sediment_mass)
                # initial_sediment_mass_wet = 0  # todo: need a better way such that we do not do this
                initial_sediment_mass_wet = np.NaN
                return initial_sediment_mass_dry, initial_sediment_mass_wet, sediment_condition
            elif 'wet' in sediment_condition:
                initial_sediment_mass_dry = np.NaN
                initial_sediment_mass_wet = float(initial_sediment_mass)
                return initial_sediment_mass_dry, initial_sediment_mass_wet, sediment_condition
            else:
                return np.NaN, np.NaN, ''


    # fetch microbial biomass: cells count in water (cells/mL water)
    def fetch_biomass_cells_count_water(self):
        try:
            raw = self.info.get_biomass_ws().get_value().split(';')[0]
            unit = self.info.get_biomass_ws().get_unit()
        except:
            return np.NaN, ''
        else:
            return raw, unit

    # fetch microbial biomass: cells count in sediment (cells/g sediment)
    def fetch_biomass_cells_count_sediment(self):
        try:
            raw = self.info.get_biomass_ws().get_value().split(';')[1]
            unit = self.info.get_biomass_ws().get_unit()
        except:
            return np.NaN, ''
        else:
            return raw, unit

    # fetch microbial biomass in sediment (mg C/g sediment)
    def fetch_biomass_sediment(self):
        try:
            raw = self.info.get_biomass_ws().get_value().split(';')[2]
            unit = self.info.get_biomass_ws().get_unit()
        except:
            return np.NaN, np.NaN, ''
        else:
            if '-' in raw:
                value, comment = self.is_censored(raw)  # comment contains about the replaced sign, if any
                value = value.split(' - ')
                if value[0] == '2023-01-1700:00:00':  # todo: remove once database is fixed
                    value[0] = 0  # temporary value, until Atorvastatin's biomass_start value is verified form DAR
                return float(value[0]), float(value[1]), unit
            elif raw == 'NA' or raw == '':
                return np.NaN, np.NaN, ''
            else:
                return np.NaN, np.NaN, ''

    # Oxygen content in water layer
    def fetch_oxygen_content_water(self):
        try:
            raw = self.info.get_oxygencontent().get_value()
            unit = self.info.get_oxygencontent().get_unit()
        except:
            return np.NaN, np.NaN, ''
        else:
            value = raw.split(';')[0]
            if ' - ' in value:
                start_value = value.split(' - ')[0]
                end_value = value.split(' - ')[1]
                if ',' in start_value or ',' in end_value:
                    start_value, var = self.is_censored(start_value)
                    end_value, var2 = self.is_censored(end_value)
                    if '(' or ')' in start_value or '(' or ')' in end_value:
                        start_value = start_value.replace('(', '')
                        start_value =start_value.replace(')', '')
                        end_value = end_value.replace('(', '')
                        end_value = end_value.replace(')', '')
                        if '-' in start_value:
                            start_value = start_value.split('-')[0]
                        if '-' in end_value:
                            end_value = end_value.split('-')[1]
                    return float(start_value), float(end_value), unit
                elif '-' in start_value or '-' in end_value:
                    start_value = start_value.replace('-', '.')
                    end_value = end_value.replace('-', '.')
                    return float(start_value), float(end_value), unit
                else:
                    return float(start_value), float(end_value), unit
            elif raw == 'NA' or raw == '':
                return np.NaN, np.NaN, ''
            else:
                return np.NaN, np.NaN, ''

    # Oxygen content in sediment
    def fetch_oxygen_content_sediment(self):
        try:
            raw = self.info.get_oxygencontent().get_value()
            unit = self.info.get_oxygencontent().get_unit()
        except:
            return np.NaN, np.NaN, ''
        else:
            value = raw.split(';')[1]
            if '-' in value:
                start_value = value.split('-')[0]
                end_value = value.split('-')[1]
                return float(start_value), float(end_value), unit
            elif raw == 'NA' or raw == '':
                return np.NaN, np.NaN, ''
            else:
                return np.NaN, np.NaN, ''

    # water-sediment functions, to be rewritten
    # fetch acidity function returns pH values in water and sediment phase, and the method [KCl, CaCl2 or H2O etc.]
    def fetch_acidity_method_sediment(self):
        try:
            raw_pH = self.info.get_acidity_ws().get_value()
        except:
            return ''
        else:
            if ';' in raw_pH:
                temp_a_method = raw_pH.split(';')[2]  # temp variable for acidity method
                if 'CaCl' in temp_a_method:
                    a_method = 'CaCl2'
                elif 'KCl' in temp_a_method:
                    a_method = 'KCl'
                elif 'H2O' or 'water' or 'WATER' in temp_a_method:  # check which other methods, and add elif as needed
                    a_method = 'H2O'
                else:  # when no string available
                    a_method = ''
                return a_method

    def fetch_acidity_water_phase(self):
        try:
            raw_pH = self.info.get_acidity_ws().get_value()
        except:
            return np.NaN
        else:
            if ';' in raw_pH:
                if ' - ' in raw_pH.split(';')[0]:
                    a_water = raw_pH.split(';')[0]
                    if ',' in a_water:
                        a_water, var = self.is_censored(a_water)
                    a_water = self.range_to_average(a_water)  # acidity in water
                else:
                    a_water = float(raw_pH.split(';')[0])
            elif '-' in raw_pH:  # if range, get mean value
                a_water = self.range_to_average(raw_pH)
            else:
                a_water = float(raw_pH)
            return np.round(a_water, 1)

    def fetch_acidity_sediment_phase(self):
        try:
            raw_pH = self.info.get_acidity_ws().get_value()
        except:
            return np.NaN
        else:
            if ';' in raw_pH:
                a_sediment = raw_pH.split(';')[1]
                if ' - ' in a_sediment:
                    a_sediment, var1 = self.is_censored(a_sediment)  # to check and replace ',' by '.'
                    if '-5' in a_sediment:  # 2 values are '6.5 - -5.15' # todo: fix this on website, as pH cannot be negative
                        a_sediment = a_sediment.replace('-5', '5')
                    a_sediment = self.range_to_average(a_sediment)  # acidity in sediment
                else:
                    a_sediment = float(raw_pH.split(';')[1])
            elif ' - ' in raw_pH:  # if range, get mean value
                a_sediment = self.range_to_average(raw_pH)
            else:
                a_sediment = float(raw_pH)
            return np.round(a_sediment, 1)

    def fetch_initial_volume_water(self):
        try:
            raw_value = self.info.get_initialvolumewater().get_value()
        except:
            return np.NaN
        else:
            return self.get_float_or_nan(raw_value)

    # fetch redox potentials of surface water
    def fetch_redox_potential_water(self):
        try:
            raw = self.info.get_redoxpotential().get_value().split(';')[0]
            raw, var = self.is_censored(raw)
        except:
            return np.NaN, np.NaN  # , np.NaN
        else:
            if ' - ' in raw:
                redox_start = raw.split(' - ')[0]
                redox_end = raw.split(' - ')[1]
                if ',' in redox_start or ',' in redox_end:
                    redox_start = redox_start.replace(',', '.')
                    redox_end = redox_end.replace(',', '.')
                return redox_start, redox_end
            else:
                return np.NaN, np.NaN

    # fetch redox potential of sediment phase
    def fetch_redox_potential_sediment(self):
        try:
            raw = self.info.get_redoxpotential().get_value().split(';')[1]
            raw, var = self.is_censored(raw)
        except:
            return np.NaN, np.NaN
        else:
            if ' - ' in raw:
                redox_start = raw.split(' - ')[0]
                redox_end = raw.split(' - ')[1]
                return redox_start, redox_end
            else:
                return raw, raw

    # fetch sample porosity
    def fetch_sample_porosity(self):
        try:
            sample_porosity = self.info.get_sedimentporosity().get_value()
        except:
            return np.NaN
        else:
            return self.get_float_or_nan(sample_porosity)

    # fetch column height for water-sediment data (cm)
    def fetch_column_height(self):
        try:
            raw_column_height = self.info.get_columnheight().get_value()
        except:
            return np.NaN, np.NaN
        else:
            column_height_sediment = self.get_float_or_nan(raw_column_height.split(';')[0])
            column_height_water = self.get_float_or_nan(raw_column_height.split(';')[1])
            return column_height_water, column_height_sediment

    # Organic carbon in water layer: Total organic carbon (TOC) and its unit
    def fetch_total_organic_carbon(self):
        try:
            raw = self.info.get_organiccarbonwater().get_value().split(';')[0]
            unit = self.info.get_organiccarbonwater().get_unit()
        except:
            return np.NaN, np.NaN, ''
        else:
            raw, var = self.is_censored(raw)
            if ' - ' in raw:
                toc1 = raw.split(' - ')[0]
                toc2 = raw.split(' - ')[1]
                if ',' in toc1 or ',' in toc2:
                    toc1, var1 = self.is_censored(toc1)
                    toc2, var2 = self.is_censored(toc2)
                if '-' in toc1 or '-' in toc2:
                    toc1 = toc1.replace('-', '.')
                    toc2 = toc2.replace('-', '.')
                return float(toc1), float(toc2), unit
            else:
                return np.NaN, np.NaN, ''

    # Organic carbon in water layer: Dissolved organic carbon (DOC) and its unit
    def fetch_dissolved_organic_carbon(self):
        try:
            raw = self.info.get_organiccarbonwater().get_value().split(';')[1]
            unit = self.info.get_organiccarbonwater().get_unit()
        except:
            return np.NaN,np.NaN, ''
        else:
            raw, var = self.is_censored(raw)
            if ' - ' in raw:
                doc1 = raw.split(' - ')[0]
                doc2 = raw.split(' - ')[1]
                return float(doc1), float(doc2), unit
            else:
                return np.NaN, np.NaN, ''

    # Organic content in sediment - organic carbon [OC] and organic matter [OM] values in %
    def fetch_organic_carbon_sediment(self):
        try:
            raw = self.info.get_organiccontent().get_value().split(';')[0]
        except:
            return np.NaN, np.NaN
        else:
            raw, var = self.is_censored(raw)
            if ' - ' in raw:
                oc1 = raw.split(' - ')[0]
                oc2 = raw.split(' - ')[1]
                oc1, var1 = self.is_censored(oc1)
                oc2, var2 = self.is_censored(oc2)
                return float(oc1), float(oc2)
            else:
                return np.NaN, np.NaN

    # High OC or Low OC scenario information
    def oc_type(self):
        try:
            oc_type = self.des
        except:
            return ''
        else:
            if 'high organic carbon content scenario' in oc_type:
                oc_type_is = 'high OC'
                return oc_type_is
            elif 'low organic carbon content scenario' in oc_type:
                oc_type_is = 'low OC'
                return oc_type_is
            else:
                return ''

    def fetch_organic_matter_sediment(self):
        try:
            raw = self.info.get_organiccontent().get_value().split(';')[1]
        except:
            return np.NaN, np.NaN
        else:
            raw, var = self.is_censored(raw)
            if ' - ' in raw:
                om1 = raw.split(' - ')[0]
                om2 = raw.split(' - ')[1]
                return float(om1), float(om2)
            else:
                return np.NaN, np.NaN