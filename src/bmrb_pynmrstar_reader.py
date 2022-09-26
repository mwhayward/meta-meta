import pynmrstar
import pandas as pd
import os
from pathlib import Path
import sqlite3
from sqlite3 import Error


class BMRB_Reader:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.conn = self.create_connection()
        self.tables = {'metabolites': {'metabolite_id': [],
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
                                   'units': [],
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

    def create_connection(self):
        # Creates a single connection to the database
        conn = None
        try:
            conn = sqlite3.connect(str(self.directory.joinpath('ccpn_metabolites_bmrb.db')))
        except Error as e:
            print(e)
        return conn

    def get_loop_tables(self, entry, category):
        # takes the pynmrstar entry object and the target category
        # returns a list of pandas dataframes
        loops = entry.get_loops_by_category(category)
        tables = [pd.DataFrame(loop.data, columns=loop.tags)for loop in loops]
        if len(tables)==0:
            pass
            #print(f'no loop table for category: {category}')
        return tables

    def get_saveframe_tags(self, entry, category):
        # takes the pynmrstar entry object and the target category
        # returns a list of pandas dataframes
        saveframes = entry.get_saveframes_by_category(category)
        saveframe_tags = [pd.DataFrame(saveframe.tags, columns=['tag', 'value']) for saveframe in saveframes]
        return saveframe_tags

    def create_database(self):
        metabolites = pd.DataFrame(self.tables['metabolites'])
        metabolites.to_sql('metabolites', self.conn, if_exists='replace', index=False)
        samples = pd.DataFrame(self.tables['samples'])
        samples.to_sql('samples', self.conn, if_exists='replace', index=False)
        spectra = pd.DataFrame(self.tables['spectra'])
        spectra.to_sql('spectra', self.conn, if_exists='replace', index=False)
        peaks = pd.DataFrame(self.tables['peaks'])
        peaks.to_sql('peaks', self.conn, if_exists='replace', index=False)

    def run(self):
        sample_added = False
        spectrum_added = False
        # put all files from the target directory into a list
        target_dir = self.directory.joinpath('BMRB_files/bmrb_nmr_spectra')
        files = os.listdir(target_dir)
        # set up counts for the reader
        metabolite_count = 1
        sample_count = 1
        spectrum_count = 1
        fail_count = 0
        # iterate through each file
        for file in files:
            filepath = target_dir.joinpath(file)
            metabolite_id = f'SU:{metabolite_count}'
            entry = pynmrstar.Entry.from_file(str(filepath))
            entry_number = entry.get_tag('Entry.ID')[0]
            if entry_number.startswith('bmse01'):
                continue
            #print(f'{metabolite_count} out of {len(files)}, ({entry_number})')
            try:
                name = entry.get_tag('Entry.Title')[0]
                # check if metabolite is already in the list
                if not name in self.tables['metabolites']['name']:
                    # populate metabolite entry if the entry is not already there
                    self.tables['metabolites']['metabolite_id'].append(metabolite_id)
                    self.tables['metabolites']['accession'].append(entry_number)
                    self.tables['metabolites']['name'].append(name)
                    self.tables['metabolites']['description'].append(None)
                    chemical_formula = entry.get_tag('Chem_comp.Formula')[0]
                    self.tables['metabolites']['chemical_formula'].append(chemical_formula)
                    molecular_weight = entry.get_tag('Chem_comp.Formula_weight')[0]
                    self.tables['metabolites']['molecular_weight'].append(molecular_weight)
                    smiles_tables = self.get_loop_tables(entry, 'Chem_comp_SMILES')
                    if len(smiles_tables)==0:
                        smiles_table = self.get_loop_tables(entry, 'Chem_comp_descriptor')[0]
                        smiles = smiles_table.loc[smiles_table['Type'] == 'SMILES', 'Descriptor'].iloc[0]
                    else:
                        smiles = smiles_tables[0].loc[smiles_tables[0]['Type'] == 'canonical', 'String'].iloc[0]
                    self.tables['metabolites']['smiles'].append(smiles)
                    inchi = entry.get_tag('Chem_comp.InChI_code')[0]
                    self.tables['metabolites']['inChi'].append(inchi)
                    metabolite_count += 1
                else:
                    # Use available entry if there is already one there
                    id_index = self.tables['metabolites']['name'].index(name)
                    metabolite_id = self.tables['metabolites']['metabolite_id'][id_index]
                # find the dimension details for all spectral peak lists
                dimension_tables = self.get_loop_tables(entry, 'Spectral_dim')
                # only use the peak lists that are 1D 1H
                peaklist_ids = []
                for table in [table for table in dimension_tables if len(table) == 1 and table.loc[0, 'Atom_type'] == 'H']:
                    peaklist_id = table.loc[table['Atom_type'] == 'H', 'Spectral_peak_list_ID'].iloc[0]
                    peaklist_ids.append(peaklist_id)
                tagtables = self.get_saveframe_tags(entry, 'spectral_peak_list')
                links = {}
                for tagtable in [tagtable for tagtable in tagtables if tagtable.loc[tagtable['tag'] == 'ID', 'value'].iloc[0] in peaklist_ids]:
                    peaklist_id_linkname = tagtable.loc[tagtable['tag'] == 'ID', 'value'].iloc[0]
                    experiment_id = tagtable.loc[tagtable['tag'] == 'Experiment_ID', 'value'].iloc[0]
                    sample_id = tagtable.loc[tagtable['tag'] == 'Sample_ID', 'value'].iloc[0]

                    links[peaklist_id_linkname] = {'experiment_id': experiment_id,
                                                   'sample_id': sample_id}
                sample_tables = self.get_loop_tables(entry, 'Sample_component')
                sample_condition_tables = self.get_loop_tables(entry, 'Sample_condition_variable')
                if len(sample_condition_tables) > 1:
                    print(f'more than one sample condition in entry {entry_number}')
                experiment_tables = self.get_loop_tables(entry, 'Experiment')
                if len(experiment_tables) > 1:
                    print(f'more than one experiment conditions table in entry {entry_number}')
                spectrometer_tags_tables = self.get_saveframe_tags(entry, 'NMR_spectrometer')
                chem_shift_tables = [table for table in self.get_loop_tables(entry, 'Spectral_transition_char') if
                                     table.loc[0, 'Spectral_peak_list_ID'] in peaklist_ids]
                intensity_tables = [table for table in self.get_loop_tables(entry, 'Spectral_transition_general_char') if
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
                                    self.tables['peaks']['peak_id'].append(peak_id)
                                    self.tables['peaks']['spectrum_id'].append(spectrum_id)
                                    self.tables['peaks']['multiplet_id'].append(None)
                                    peak_shift = row['Chem_shift_val']
                                    self.tables['peaks']['shift'].append(peak_shift)
                                    peak_intensity = intensity_tables[table_num].loc[index, 'Intensity_val']
                                    self.tables['peaks']['intensity'].append(peak_intensity)
                                    self.tables['peaks']['width'].append(None)
                                if sample_added is False:
                                    self.tables['samples']['sample_id'].append(sample_id)
                                    self.tables['samples']['metabolite_id'].append(metabolite_id)
                                    self.tables['samples']['pH'].append(ph)
                                    self.tables['samples']['temperature'].append(temperture)
                                    self.tables['samples']['amount'].append(amount)
                                    self.tables['samples']['units'].append(units)
                                    self.tables['samples']['reference'].append(reference)
                                    sample_count += 1
                                    sample_added = True
                                if spectrum_added is False:
                                    self.tables['spectra']['spectrum_id'].append(spectrum_id)
                                    self.tables['spectra']['sample_id'].append(sample_id)
                                    self.tables['spectra']['frequency'].append(frequency)
                                    spectrum_count += 1
                                    spectrum_added = True
            except:
                print(f'Unable to parse entry {entry_number}')
                fail_count +=1
        print(f'fail count = {fail_count}')


if __name__ == "__main__":
    directory = '/home/mh491/Database'
    reader = BMRB_Reader(directory)
    reader.run()
    reader.create_database()
