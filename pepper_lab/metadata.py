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
        D = {'compound_id': [], 'compound_name': [], 'smiles': [],  # compound
             'scenario_id': [], 'study_name': [], 'study_description': [], # study/scenario
             'halflife_raw': [], 'halflife_unit': [], 'halflife_model': [], 'halflife_comment': [], # DT50
             'spike_compound': [],
             # additional information
             'acidity': [], 'acidity_unit': [],
             'temperature': [], 'temperature_unit': [],
             'CEC': [],
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
            "compound_id": [], "compound_name": [], "smiles": [], # compound
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
        D = {'compound_id': [], 'compound_name': [], 'smiles': [], # compound information
             'scenario_id': [], 'study_name': [], 'study_description': [], # study description
             # DT50
             'DT50_water': [], 'DT50_sediment': [], 'DT50_total_system': [],
             'DT50_water_comment': [], 'DT50_sediment_comment': [], 'DT50_total_system_comment': [],
             'halflife_model': [], 'halflife_comment': [], 'halflife_fit': [],
             'spike_compound': [],
             # additional information
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
             'biomass_sediment_start': [], 'biomass_sediment_end': [],  # range of biomass at start and end in sediment
             'biomass': [],  # avg biomass (average of biomass_sediment_start and biomass_sediment_end)
             'biomass_sediment_unit': [],
             'temperature': [], 'temperature_unit':[],
             'sample_location': [],
             'sediment_porosity': [],
             'initial_sediment_mass_dry': [], 'initial_sediment_mass_wet': [],
             'sediment_condition': [],
             'initial_volume_water': [],
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

    def get_scenario_information(self, D, scenario, compound, data_type, spike_smiles, description):
        if data_type == 'soil':
            D = self.get_soil_scenario_information(D, scenario, compound, spike_smiles, description)
        elif data_type == 'sediment':
            D = self.get_sediment_scenario_information(D, scenario, compound, spike_smiles, description)
        elif data_type == 'sludge':
            D = self.get_sludge_scenario_information(D, scenario, compound, description)
        else:
            raise NotImplementedError
        return D

    def get_sludge_scenario_information(self, D, scenario, compound, description):
        if D == {}:
            D = self.initiate_sludge_dictionary()
        # Compound informatin
        D['compound_id'].append(compound.get_id())
        D['compound_name'].append(compound.get_name())
        D['smiles'].append(compound.get_smiles())
        # Scenario/study information
        D['scenario_id'].append(scenario.get_id())
        D['study_name'].append(scenario.get_name().split(' - ')[0])
        D['study_description'].append(description)
        D['acidity'].append(self.fetch_mean_value("acidity", 'lowPh', 'highPh'))
        D['acidity_unit'].append(self.fetch_normal_value("acidity", "unit", str))
        D['addition_of_nutrients'].append(self.fetch_normal_value("additionofnutrients", "additionofnutrients", str))
        D['biological_treatment_technology'].append(self.fetch_normal_value("biologicaltreatmenttechnology", "biologicaltreatmenttechnology", str))
        D['bioreactor_type'].append(self.fetch_normal_value("bioreactor", "bioreactortype", str))
        D['bioreactor_value'].append(self.fetch_normal_value("bioreactor", "bioreactorsize", float))
        D['bioreactor_value_unit'].append(self.fetch_normal_value("bioreactor", "unit", str))
        D['halflife_raw'].append(self.fetch_mean_value("halflife", "lower", "upper"))
        D['halflife_unit'].append(self.fetch_normal_value("halflife", "unit", str))
        D['halflife_model_TF'].append(self.fetch_normal_value("halflife", "model", str))
        D['halflife_comment'].append(self.fetch_normal_value("halflife", "comment", str))
        D['inoculum_source'].append(self.fetch_normal_value("inoculumsource", "inoculumsource", str))
        D['location'].append(self.fetch_normal_value("location", "location", str))
        D['nitrogen_content_type'].append(self.fetch_normal_value("nitrogencontent", "nitrogencontentType", str))
        D['nitrogen_content_influent'].append(self.fetch_normal_value("nitrogencontent", "nitrogencontentInfluent", str))
        D['original_sludge_amount'].append(self.fetch_normal_value('originalsludgeamount', 'originalsludgeamount', float))
        D['original_sludge_amount_unit'].append(self.fetch_normal_value('originalsludgeamount', 'unit', str))
        D['oxygen_demand_type'].append(self.fetch_normal_value('oxygendemand', 'oxygendemandType', str))
        D['oxygen_demand_value'].append(self.fetch_normal_value('oxygendemand', 'oxygendemandInfluent', float))
        D['oxygen_uptake_rate_unit'].append(self.fetch_normal_value("oxygenuptakerate", "unit", str))
        D['oxygen_uptake_rate'].append(self.fetch_mean_value('oxygenuptakerate', 'oxygenuptakerateStart', 'oxygenuptakerateEnd'))
        D['phosphorus_content'].append(self.fetch_normal_value("phosphoruscontent", "phosphoruscontentInfluent", float))
        D['purpose_of_wwtp'].append(self.fetch_normal_value("purposeofwwtp", "purposeofwwtp", str))
        D['rateconstant'].append(self.fetch_mean_value("rateconstant", "rateconstantlower", "rateconstantupper"))
        D['rateconstant_unit'].append(self.fetch_normal_value("rateconstant", "unit", str))
        D['halflife_model'].append(self.fetch_normal_value("rateconstant", 'rateconstantorder', str))
        D['rateconstant_comment'].append(self.fetch_normal_value("rateconstant", 'rateconstantcomment', str))
        D['redox'].append(self.fetch_normal_value("redox", "redoxType", str))
        D['sludge_retention_time_type'].append(self.fetch_normal_value("sludgeretentiontime", "sludgeretentiontimeType", str))
        D['sludge_retention_time'].append(self.fetch_normal_value("sludgeretentiontime", "sludgeretentiontime", float))
        D['sludge_retention_time_unit'].append(self.fetch_normal_value("sludgeretentiontime", "unit", str))
        D['source_of_liquid_matrix'].append(self.fetch_normal_value("sourceofliquidmatrix", "sourceofliquidmatrix", str))
        D['temperature'].append(self.fetch_mean_value("temperature", "temperatureMin", "temperatureMax"))
        D['temperature_unit'].append(self.fetch_normal_value("temperature", "unit", str))
        D['total_suspended_solids_concentration_start'].append(self.fetch_normal_value("tts", "ttsStart", float))
        D['total_suspended_solids_concentration_end'].append(self.fetch_normal_value("tts", "ttsEnd", float))
        D['total_suspended_solids_concentration_unit'].append(self.fetch_normal_value("tts", "unit", str))
        D['type_of_addition'].append(self.fetch_normal_value("typeofaddition", "typeofaddition", str))
        D['type_of_aeration'].append(self.fetch_normal_value("typeofaeration", "typeofaeration", str))
        return D

    def get_soil_scenario_information(self, D, scenario, compound, spike_smiles, description):
        # compound info
        if D == {}:
            D = self.initiate_soil_dictionary()
        D['compound_id'].append(compound.get_id())
        D['compound_name'].append(compound.get_name())
        D['smiles'].append(compound.get_smiles())
        D['spike_compound'].append(spike_smiles)
        # study
        D['scenario_id'].append(scenario.get_id())
        D['study_name'].append(scenario.get_name().split(' - ')[0])
        D['study_description'].append(description)

        # add halflife details
        D['halflife_raw'].append(self.fetch_mean_value("halflife", "lower", "upper")) # renamed from 'reported_DT50'
        D['halflife_unit'].append(self.fetch_normal_value("halflife", "unit", str)) # added
        D['halflife_model'].append(self.fetch_normal_value("halflife", "model", str))
        D['halflife_comment'].append(self.fetch_normal_value("halflife", "comment", str))

        # fetch additional information
        D['acidity'].append(self.fetch_mean_value("acidity", 'lowPh', 'highPh'))
        D['acidity_unit'].append(self.fetch_normal_value("acidity", "unit", str)) # added
        D['temperature'].append(self.fetch_mean_value("temperature", "temperatureMin", "temperatureMax"))
        D['temperature_unit'].append(self.fetch_normal_value("temperature", "unit", str))

        D['CEC'].append(self.fetch_cec())  # cation exchange capacity
        D['OC'].append(self.fetch_organic_content())  # organic content as organic carbon (oc)
        start, end = self.fetch_biomass()
        D['biomass_start'].append(start)
        D['biomass_end'].append(end)
        D['biomass'].append(np.round(np.average([start, end]), 2))
        wst_value, wst_type = self.fetch_wst()  # water storage capacity,
        D['wst_value'].append(wst_value)
        D['wst_type'].append(wst_type)
        hum, hum_cond = self.fetch_humidity()
        D['humidity'].append(hum)
        D['humidity_conditions'].append(hum_cond)
        D['soil_texture'].append(self.fetch_soiltexture1())
        _sand, _silt, _clay = self.fetch_soiltexture2()
        D['sand'].append(_sand)
        D['silt'].append(_silt)
        D['clay'].append(_clay)
        return D

    def get_sediment_scenario_information(self, D, scenario, compound, spike_smiles, description):
        # compound info
        if D == {}:
            D = self.initiate_sediment_dictionary()
        D['compound_id'].append(compound.get_id())
        D['compound_name'].append(compound.get_name())
        D['smiles'].append(compound.get_smiles())
        D['scenario_id'].append(scenario.get_id())
        D['study_description'].append(description)
        # fetch halflife details - total system, water, sediment, model_type, comment, fit
        dt50_total, comment1 = self.fetch_halflife_total_system_value()
        D['DT50_total_system'].append(self.range_to_average(dt50_total))
        D['DT50_total_system_comment'].append(comment1)
        dt50_water, comment2 = self.fetch_halflife_water_value()
        D['DT50_water'].append(self.range_to_average(dt50_water))
        D['DT50_water_comment'].append(comment2)
        dt50_sediment, comment3 = self.fetch_halflife_sediment_value()
        D['DT50_sediment'].append(self.range_to_average(dt50_sediment))
        D['DT50_sediment_comment'].append(comment3)
        D['halflife_model'].append(self.fetch_halflife_ws_model())
        D['halflife_comment'].append(self.fetch_halflife_ws_comment())
        D['halflife_fit'].append(self.fetch_halflife_ws_fit())
        # Fetch other data points
        D['study_name'].append(scenario.get_name().split(' - ')[0])
        D['spike_compound'].append(spike_smiles)
        #  fetch pH values for surface water and sediment, and the method used for measuring pH in sediment
        D['temperature'].append(self.fetch_mean_value("temperature", "temperatureMin", "temperatureMax"))
        D['temperature_unit'].append(self.fetch_normal_value("temperature", "unit", str))
        D['acidity_water'].append(self.fetch_acidity_water_phase())  # pH surface water
        D['acidity_sediment'].append(self.fetch_acidity_sediment_phase())  # pH sediment
        D['acidity_method'].append(self.fetch_acidity_method_sediment())  # pH method in sediment
        D['bulk_density'].append(self.fetch_bulk_density())  # bulk density
        D['bulk_density_unit'].append(self.fetch_bulk_density_unit())  # unit of bulk density
        D['CEC'].append(self.fetch_cec())  # cation exchange capacity
        # column height for water and sediment phases, respectively
        column_height_w, column_height_s = self.fetch_column_height()
        D['column_height_water'].append(column_height_w)
        D['column_height_sediment'].append(column_height_s)
        # initial sediment mass (dry/wet)
        initial_sediment_mass_d, initial_sediment_mass_w, sediment_condition = self.fetch_initial_sediment_mass()
        D['initial_sediment_mass_dry'].append(initial_sediment_mass_d)
        D['initial_sediment_mass_wet'].append(initial_sediment_mass_w)
        D['sediment_condition'].append(sediment_condition)  # dry or wet
        # initial volume of water
        D['initial_volume_water'].append(self.fetch_initial_volume_water())
        # System information: High OC or Low OC system
        D['OC_type'].append(self.oc_type())  # high oc / low oc vales
        # Organic carbon in water layer - Total Organic Carbon [TOC] values
        toc1, toc2, toc_unit = self.fetch_total_organic_carbon()
        D['TOC_1'].append(toc1)
        D['TOC_2'].append(toc2)
        D['TOC'].append(np.round(np.average([toc1, toc2]), 2))
        D['TOC_unit'].append(toc_unit)
        # Organic carbon in water layer - Dissolved Organic Carbon [DOC] values
        doc1, doc2, doc_unit = self.fetch_dissolved_organic_carbon()
        D['DOC_1'].append(doc1)
        D['DOC_2'].append(doc2)
        D['DOC'].append(np.round(np.average([doc1, doc2]), 2))
        D['DOC_unit'].append(doc_unit)
        # Organic content in sediment organic carbon [OC] and organic matter [OM] values
        oc1, oc2 = self.fetch_organic_carbon_sediment()
        D['OC_1'].append(oc1)
        D['OC_2'].append(oc2)
        D['OC'].append(np.round(np.average([oc1, oc2]), 2))
        om1, om2 = self.fetch_organic_matter_sediment()
        D['OM_1'].append(om1)
        D['OM_2'].append(om2)
        D['OM'].append(np.round(np.average([om1, om2]), 2))
        # Oxygen content in water layer
        start_oxygen_water, end_oxygen_water, oxygen_water_unit = self.fetch_oxygen_content_water()
        D['oxygen_content_water_start'].append(start_oxygen_water)
        D['oxygen_content_water_end'].append(end_oxygen_water)
        D['oxygen_content_water'].append(np.round(np.average([start_oxygen_water, end_oxygen_water]), 2))  # avg
        D['oxygen_content_water_unit'].append(oxygen_water_unit)
        # Oxygen content in sediment layer
        oxygen_sediment_start, oxygen_sediment_end, oxygen_sediment_unit = self.fetch_oxygen_content_sediment()
        D['oxygen_content_sediment_start'].append(oxygen_sediment_start)
        D['oxygen_content_sediment_end'].append(oxygen_sediment_end)
        D['oxygen_content_sediment'].append(np.round(np.average([oxygen_sediment_start, oxygen_sediment_end]), 2))
        D['oxygen_content_sediment_unit'].append(oxygen_sediment_unit)
        # Microbial biomass in water
        cells_water_count, cells_water_count_unit = self.fetch_biomass_cells_count_water()
        D['biomass_cells_count_water'].append(cells_water_count)
        D['biomass_cells_count_water_unit'].append(cells_water_count_unit)
        # Microbial biomass in sediment
        cells_count_sediment, cells_count_sediment_unit = self.fetch_biomass_cells_count_sediment()
        D['biomass_cells_count_sediment'].append(cells_count_sediment)
        D['biomass_cells_count_sediment_unit'].append(cells_count_sediment_unit)
        start_sediment, end_sediment, biomass_unit = self.fetch_biomass_sediment()
        D['biomass_sediment_start'].append(start_sediment)
        D['biomass_sediment_end'].append(end_sediment)
        D['biomass'].append(np.round(np.average([start_sediment, end_sediment]), 2))
        D['biomass_sediment_unit'].append(biomass_unit)
        # Redox potentials
        redox_w1, redox_w2 = self.fetch_redox_potential_water()  # Redox potential of water at start and end, respectively
        D['redox_water_start'].append(redox_w1)
        D['redox_water_end'].append(redox_w2)
        redox_s1, redox_s2 = self.fetch_redox_potential_sediment()  # Redox potential of sediment at start and end, respectively
        D['redox_sediment_start'].append(redox_s1)
        D['redox_sediment_end'].append(redox_s2)
        # sample location
        D['sample_location'].append(self.fetch_sample_location())
        D['sediment_porosity'].append(self.fetch_sample_porosity())  # sediment porosity
        D['soil_texture'].append(self.fetch_soiltexture1())
        _sand, _silt, _clay = self.fetch_soiltexture2()
        D['sand'].append(_sand)
        D['silt'].append(_silt)
        D['clay'].append(_clay)
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