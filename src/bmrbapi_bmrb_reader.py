import pandas as pd
from pathlib import Path
import sqlite3
from sqlite3 import Error
import requests
import json


def get_bmrb_entry(entry_number):
    # takes entry number
    # returns a json of the whole entry
    data = requests.get(f"http://api.bmrb.io/v2/entry/{entry_number}", headers={"Application": "CCPN_Analysis_Metablomics V3"})
    data_out = json.loads(data.text)
    return data_out


def get_bmrb_tag(entry_number, tag):
    # takes entry number and the target tag
    # returns a string result
    data = requests.get(f"http://api.bmrb.io/v2/entry/{entry_number}?tag={tag}", headers={"Application": "CCPN_Analysis_Metablomics V3"})
    data_out = json.loads(data.text)[entry_number][tag][0]
    return data_out


def get_bmrb_save(entry_number, category):
    # takes entry number and saveframe category
    # returns a list of saveframes in json format
    data = requests.get(f"http://api.bmrb.io/v2/entry/{entry_number}?saveframe_category={category}", headers={"Application": "CCPN_Analysis_Metablomics V3"})
    json_data = json.loads(data.text)[entry_number][category]
    saveframes = [saveframe for saveframe in json_data]
    return saveframes


def get_bmrb_save_tags(entry_number, category):
    # takes entry number and saveframe category
    # returns a list of saveframe tags in pandas dataframes
    data = requests.get(f"http://api.bmrb.io/v2/entry/{entry_number}?saveframe_category={category}", headers={"Application": "CCPN_Analysis_Metablomics V3"})
    json_data = json.loads(data.text)[entry_number][category]
    saveframe_tags = [pd.DataFrame(saveframe['tags'], columns=['tag', 'value']) for saveframe in json_data]
    return saveframe_tags


def get_bmrb_loop(entry_number, tag):
    # takes entry number and loop prefix (without underscore)
    # returns a list of pandas dataframes for each loop recovered
    data = requests.get(f"http://api.bmrb.io/v2/entry/{entry_number}?loop={tag}", headers={"Application": "CCPN_Analysis_Metablomics V3"})
    json_data = json.loads(data.text)
    loops = []
    for loop in json_data[entry_number][tag]:
        titles = loop['tags']
        table_data = loop['data']
        data_out = pd.DataFrame(table_data, columns=titles)
        loops.append(data_out)
    return loops


def check_database(name, conn):
    # check the database if the entry already exists
    # currently only works with name
    return False
    query = f"SELECT * FROM metabolites WHERE name is '{name}';"
    data = pd.read_sql(query, conn)
    return len(data) > 0


def get_database_entry_by_name(name, conn):
    query = f"SELECT metabolite_id FROM metabolites WHERE name = '{name}';"
    data = pd.read_sql(query, conn)
    metabolite_id = data.loc[0, 'metabolite_id']
    return metabolite_id


def parse_chemical_shift_saveframe(entry_number):
    # takes the entry number of the bmrb file and scans through the save frames for the 1D data
    # returns chemical shift, sample and spectrum values or just puts them straight into the tables
    get_bmrb_loop(entry_number, 'Spectral_dim')
    pass


def create_database(tables):
    metabolites = pd.DataFrame(tables['metabolites'])
    metabolites.to_sql('metabolites', conn, if_exists='replace', index=False)
    samples = pd.DataFrame(tables['samples'])
    samples.to_sql('samples', conn, if_exists='replace', index=False)
    spectra = pd.DataFrame(tables['spectra'])
    spectra.to_sql('spectra', conn, if_exists='replace', index=False)
    peaks = pd.DataFrame(tables['peaks'])
    peaks.to_sql('peaks', conn, if_exists='replace', index=False)


# Creates a single connection to the database
directory = Path('/home/mh491/Database/ccpn_metabolites_bmrb.db')
conn = None
try:
    conn = sqlite3.connect(str(directory))
except Error as e:
    print(e)

# establish tables
tables = {'metabolites': {'metabolite_id': [],
                          'accession': [],
                          'name': [],
                          'description': [],
                          'chemical_formula': [],
                          'molecular_weight': [],
                          'smiles': [],
                          'inChi': []},
          'samples': {'sample_id': [],
                      'metabolite_id': [],
                      'pH': [],
                      'temperature': [],
                      'amount': [],
                      'reference': []},
          'spectra': {'spectrum_id': [],
                      'sample_id': [],
                      'frequency': []},
          'multiplets': {'multiplet_id': [],
                         'spectrum_id': [],
                         'center': [],
                         'atom_ref': [],
                         'multiplicity': []},
          'peaks': {'peak_id': [],
                    'spectrum_id': [],
                    'multiplet_id': [],
                    'shift': [],
                    'intensity': [],
                    'width': []},
          'synonyms': {'metabolite_id': [],
                       'synonym': []},
          'isin': {'metabolite_id': [],
                   'group': []},
          'ontology': {'group': [],
                       'definition': []}}

# reset the database
create_database(tables)

# gather the list of bmrb entries
entries_request = requests.get("http://api.bmrb.io/v2/list_entries?database=metabolomics", headers={"Application": "CCPN_Analysis_Metablomics V3"})
entries = json.loads(entries_request.text)

# establish counts
metabolite_count = 1
sample_count = 1
spectrum_count = 1

# parse each experimental entry
for number, entry_number in enumerate([entry for entry in entries if entry.startswith('bmse')]):
    sample_added = False
    spectrum_added = False
    metabolite_id = f'SU:{metabolite_count}'
    name = get_bmrb_tag(entry_number, 'Entry.Title')
    print(f'{metabolite_count} out of {len(entries)}, ({entry_number})')

    # check if metabolite is already in the list
    if not check_database(name, conn):
        # populate metabolite entry if the entry is not already there
        tables['metabolites']['metabolite_id'].append(metabolite_id)
        tables['metabolites']['accession'].append(entry_number)
        tables['metabolites']['name'].append(name)
        tables['metabolites']['description'].append(None)
        chemical_formula = get_bmrb_tag(entry_number, 'Chem_comp.Formula')
        tables['metabolites']['chemical_formula'].append(chemical_formula)
        molecular_weight = get_bmrb_tag(entry_number, 'Chem_comp.Formula_weight')
        tables['metabolites']['molecular_weight'].append(molecular_weight)
        smiles_table = get_bmrb_loop(entry_number, 'Chem_comp_SMILES')[0]
        smiles = smiles_table.loc[smiles_table['Type'] == 'canonical', 'String'].iloc[0]
        tables['metabolites']['smiles'].append(smiles)
        inchi = get_bmrb_tag(entry_number, 'Chem_comp.InChI_code')
        tables['metabolites']['inChi'].append(inchi)
        metabolite_count += 1
    else:
        # Use available entry if there is already one there
        metabolite_id = get_database_entry_by_name(name, conn)

    # find the dimension details for all spectral peak lists
    dimension_tables = get_bmrb_loop(entry_number, 'Spectral_dim')

    # only use the peak lists that are 1D 1H
    peaklist_ids = []
    for table in [table for table in dimension_tables if len(table) == 1 and table.loc[0, 'Atom_type'] == 'H']:
        peaklist_id = table.loc[table['Atom_type'] == 'H', 'Spectral_peak_list_ID'].iloc[0]
        peaklist_ids.append(peaklist_id)
    tagtables = get_bmrb_save_tags(entry_number, 'spectral_peak_list')
    links = {}
    for tagtable in [tagtable for tagtable in tagtables if tagtable.loc[tagtable['tag'] == 'ID', 'value'].iloc[0] in peaklist_ids]:
        peaklist_id_linkname = tagtable.loc[tagtable['tag'] == 'ID', 'value'].iloc[0]
        experiment_id = tagtable.loc[tagtable['tag'] == 'Experiment_ID', 'value'].iloc[0]
        sample_id = tagtable.loc[tagtable['tag'] == 'Sample_ID', 'value'].iloc[0]

        links[peaklist_id_linkname] = {'experiment_id': experiment_id,
                                       'sample_id': sample_id}
    sample_tables = get_bmrb_loop(entry_number, 'Sample_component')
    sample_condition_tables = get_bmrb_loop(entry_number, 'Sample_condition_variable')
    if len(sample_condition_tables) > 1:
        print(f'more than one sample condition in entry {entry_number}')
    experiment_tables = get_bmrb_loop(entry_number, 'Experiment')
    if len(experiment_tables) > 1:
        print(f'more than one experiment conditions table in entry {entry_number}')
    spectrometer_tags_tables = get_bmrb_save_tags(entry_number, 'NMR_spectrometer')
    chem_shift_tables = [table for table in get_bmrb_loop(entry_number, 'Spectral_transition_char') if
                         table.loc[0, 'Spectral_peak_list_ID'] in peaklist_ids]
    intensity_tables = [table for table in get_bmrb_loop(entry_number, 'Spectral_transition_general_char') if
                        table.loc[0, 'Spectral_peak_list_ID'] in peaklist_ids]
    for sample_table in sample_tables:
        sample_id = f'SA:{sample_count}'
        bmrb_sample_id = sample_table.loc[0, 'Sample_ID']
        amount = sample_table.loc[sample_table['Type'] == 'Solute', 'Concentration_val'].iloc[0]
        units = sample_table.loc[sample_table['Type'] == 'Solute', 'Concentration_val_units'].iloc[0]
        reference = sample_table.loc[sample_table['Type'] == 'Reference', 'Mol_common_name'].iloc[0]
        ph = sample_condition_tables[0].loc[sample_condition_tables[0]['Type'] == 'pH', 'Val'].iloc[0]
        temperture = sample_condition_tables[0].loc[sample_condition_tables[0]['Type'] == 'temperature', 'Val'].iloc[0]
        for index, row in experiment_tables[0].iterrows():
            spectrum_id = f'SP:{spectrum_count}'
            experiment_id = row['ID']
            spectrometer_id = row['NMR_spectrometer_ID']
            for spectrometer_tags_table in spectrometer_tags_tables:
                if spectrometer_tags_table.loc[spectrometer_tags_table['tag'] == 'ID', 'value'].iloc[0] == spectrometer_id:
                    frequency = spectrometer_tags_table.loc[spectrometer_tags_table['tag'] == 'Field_strength', 'value'].iloc[0]
                else:
                    frequency = None

            # obtain chemical shift and intensity data from BMRB (filtered by experiment type)
            for table_num, table in enumerate(chem_shift_tables):
                peaklist_id = table.loc[0, 'Spectral_peak_list_ID']
                if links[peaklist_id]['experiment_id'] == experiment_id and links[peaklist_id]['sample_id'] == bmrb_sample_id:
                    # get the peak data
                    for index, row in table.iterrows():
                        peak_id = f'PK:{spectrum_count}.{index+1}'
                        tables['peaks']['peak_id'].append(peak_id)
                        tables['peaks']['spectrum_id'].append(spectrum_id)
                        tables['peaks']['multiplet_id'].append(None)
                        peak_shift = row['Chem_shift_val']
                        tables['peaks']['shift'].append(peak_shift)
                        peak_intensity = intensity_tables[table_num].loc[index, 'Intensity_val']
                        tables['peaks']['intensity'].append(peak_intensity)
                        tables['peaks']['width'].append(None)
                    if sample_added is False:
                        tables['samples']['sample_id'].append(sample_id)
                        tables['samples']['metabolite_id'].append(metabolite_id)
                        tables['samples']['pH'].append(ph)
                        tables['samples']['temperature'].append(temperture)
                        tables['samples']['amount'].append(amount)
                        tables['samples']['reference'].append(reference)
                        sample_count += 1
                        sample_added = True
                    if spectrum_added is False:
                        tables['spectra']['spectrum_id'].append(spectrum_id)
                        tables['spectra']['sample_id'].append(sample_id)
                        tables['spectra']['frequency'].append(frequency)
                        spectrum_count += 1
                        spectrum_added = True
    # todo remove break when successfully parse a single entry
    if metabolite_count >9:
        break
create_database(tables)

