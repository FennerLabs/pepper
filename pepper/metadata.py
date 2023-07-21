import numpy as np
import re

class Metadata:
    '''
    This class handles additional information from enviPath scenarios
    '''
    def __init__(self, additional_info):
        self.info = additional_info

    @staticmethod
    def range_to_average(input_string):
        if type(input_string) == float: #  in case we get NaN here
            return input_string
        min = float(input_string.split(' - ')[0])
        max = float(input_string.split(' - ')[1])
        avg = np.average([min, max])
        return avg

    @staticmethod
    def initiate_dictionary():
        D = {'compound_id': [], 'compound_name': [], 'smiles': [], 'reported_DT50': [], 'scenario_id': [],
             'study_name': [],
             'halflife_model': [], 'halflife_comment': [], 'spike_compound': [], 'acidity': [], 'CEC': [], 'OC': [],
             'biomass_start': [], 'biomass_end': [], 'biomass': [], 'temperature': [], 'wst_value': [],
             'wst_type': [], 'humidity': [], 'humidity_conditions': [], 'soil_texture': [], 'sand': [], 'silt': [],
             'clay': []}
        return D

    def fetch_acidity(self):
        try:
            raw_pH = self.info.get_acidity().get_value()
        except:
            return np.NaN
        else:
            if ';' in raw_pH:
                if '-' in raw_pH.split(';')[0]:
                    pH = self.range_to_average(raw_pH.split(';')[0])
                else:
                    pH = float(raw_pH.split(';')[0])
            elif '-' in raw_pH:  # if range, get mean value
                pH = self.range_to_average(raw_pH)
            else:
                pH = float(raw_pH)
            return np.round(pH, 1)

    def fetch_cec(self):
        try:
            cec = self.info.get_cec().get_value()
        except:
            return np.NaN
        else:
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
            min = float(raw.split(';')[0])
            max = float(raw.split(';')[1])
            return np.round(np.average([min, max]), 0)

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
            values = re.findall(r'\s([\d.]+)%', raw)
            if values == []:
                return np.NaN, np.NaN, np.NaN
            return self.get_float_or_nan(values[0]), self.get_float_or_nan(values[1]), self.get_float_or_nan(
                values[2])  # sand, silt, clay

    def fetch_halflife_value(self):
        try:
            raw = self.info.get_halflife().get_value()
        except:
            return ''
        else:
            if len(raw.split(';')) < 4:
                print('Warning: incomplete half-life information - {}'.format(raw))
                val = np.NaN
            else:
                val = raw.split(';')[3]
            return val

    def fetch_halflife_model(self):
        try:
            raw = self.info.get_halflife().get_value()
        except:
            return ''
        else:
            return raw.split(';')[0]


    def fetch_halflife_comment(self):
        try:
            raw = self.info.get_halflife().get_value()
        except:
            return ''
        else:
            return raw.split(';')[2]

    @staticmethod
    def get_float_or_nan(x):
        try:
            return float(x)
        except:
            return np.NaN

    def get_scenario_information(self, D, scenario, compound, data_type, spike_smiles):
        if data_type == 'soil':
            D = self.get_soil_scenario_information(D, scenario, compound, spike_smiles)
        else:
            raise NotImplementedError
        return D

    def get_soil_scenario_information(self, D, scenario, compound, spike_smiles):
        # compound info
        if D == {}:
            D = self.initiate_dictionary()
        D['compound_id'].append(compound.get_id())
        D['compound_name'].append(compound.get_name())
        D['smiles'].append(compound.get_smiles())
        D['scenario_id'].append(scenario.get_id())
        # add halflife details
        D['reported_DT50'].append(self.range_to_average(self.fetch_halflife_value()))
        D['halflife_model'].append(self.fetch_halflife_model())
        D['halflife_comment'].append(self.fetch_halflife_comment())
        D['study_name'].append(scenario.get_name().split(' - ')[0])
        D['spike_compound'].append(spike_smiles)
        # fetch data points
        D['acidity'].append(self.fetch_acidity())
        D['CEC'].append(self.fetch_cec())  # cation exchange capacity
        D['OC'].append(self.fetch_organic_content())  # organic content as organic carbon (oc)
        start, end = self.fetch_biomass()
        D['biomass_start'].append(start)
        D['biomass_end'].append(end)
        D['biomass'].append(np.round(np.average([start, end]), 2))
        D['temperature'].append(self.fetch_temperature())
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