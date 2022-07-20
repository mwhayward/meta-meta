from pathlib import Path
import xml.etree.ElementTree as et
import nmrstarlib
import nmrpystar
import os
import sqlite3
from sqlite3 import Error

import pandas as pd


class BMRB_Reader:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.conn = self.create_connection()

    def create_connection(self):
        # Creates a single connection to the database
        conn = None
        try:
            conn = sqlite3.connect(str(self.directory.joinpath('ccpn_metabolites.db')))
        except Error as e:
            print(e)
        return conn

    def run(self):
        count = 0
        # define the dictionaries for populating
        metabolites = {'metabolite_id': [],
                       'name': [],
                       'description': [],
                       'chemical_formula': [],
                       'molecular_weight': [],
                       'smiles': [],
                       'inChi': []}
        samples = {'sample_id': [],
                   'metabolite_id': [],
                   'pH': [],
                   'amount': [],
                   'reference': []}
        spectra = {'spectrum_id': [],
                   'sample_id': [],
                   'temperature': [],
                   'frequency': []}
        multiplets = {'multiplet_id': [],
                      'spectrum_id': [],
                      'center': [],
                      'atom_ref': [],
                      'multiplicity': []}
        peaks = {'peak_id': [],
                 'spectrum_id': [],
                 'multiplet_id': [],
                 'shift': [],
                 'intensity': [],
                 'width': []}
        synonyms = {'metabolite_id': [],
                    'synonym': []}
        isin = {'metabolite_id': [],
                'group': []}
        ontology = {'group': [],
                    'definition': []}

        # dictionaries for search terms in BMSE files
        targets = {'name': ['save_assembly_1'],
                   'chemical_formula': ['Chem_comp.Formula'],
                   'smiles': ['Chem_comp_SMILES'],
                   'inChi': ['Chem_comp.InChI_code'],
                   'pH': ['save_sample_conditions_1', 'save_conditions_1', 'save_conditions_1H_DW'],
                   'amount': ['save_sample_1', 'save_sample_1H_050115_P00_02_IS_DW'],
                   'natural_source': ['save_natural_source'],
                   'peaks': ['save_spectral_peak_1H', 'save_spectral_peaks_1D_1H',
                             'save_spectral_peaks_1D_1H_set01', 'save_spectral_peaks_1D_1H_set02',
                             'save_spectral_peak_1Hp5', 'save_spectral_peak_1H_2', 'save_spectral_peak_1H_2_pH4p13',
                             'save_spectral_peak_1Hp5_pH4.13', 'save_spectral_peak_1H_pH7p5',
                             'save_spectral_peak_1H_pH6p14', 'save_spectral_peak_1H_pH8p01',
                             'save_spectral_peak_1H_2_pH4p64', 'save_spectral_peak_1Hp5_pH4.64']}
        exceptions = {'peaks': ['13C', 'DEPT', 'TOCSY', 'HSQC', 'HMBC', 'COSY']}
        target_dir = self.directory.joinpath('BMRB_files/bmrb_nmr_spectra')
        files = os.listdir(target_dir)
        for num, file in enumerate(files):
            if file.endswith('.str'):
                print(f'{num + 1} out of {len(files)}:    {file}')
                for tree in nmrstarlib.read_files(str(target_dir.joinpath(file))):
                    if 'NMR quality control of fragment libraries for screening' in tree['save_entry_information']['Entry.Title']:
                        name = tree['save_assembly_1']['Assembly.Name'].strip('\n')
                    else:
                        name = tree['save_entry_information']['Entry.Title'].strip('\n')
                    if name not in metabolites['name']:
                        metabolites['metabolite_id'].append(count)
                        count += 1
                        metabolites['name'].append(name)
                        metabolites['description'].append(None)
                        metabolites['chemical_formula'].append(self.find_save_data(tree, 'Chem_comp.Formula'))
                        metabolites['molecular_weight'].append(self.find_save_data(tree, 'Chem_comp.Formula_weight'))
                        metabolites['smiles'].append(self.find_save_data(tree, 'Chem_comp_SMILES'))
                        metabolites['inChi'].append(self.find_save_data(tree, 'Chem_comp.InChI_code'))
        metabolites_df = pd.DataFrame(metabolites)
        print(metabolites_df)

    def find_save_data(self, tree, target):
        # returns the first instance of save data that matches the given target
        # doesn't delve into loops
        for save in tree.keys():
            if not isinstance(tree[save], str):
                for entity in tree[save].keys():
                    if target == entity:
                        return tree[save][entity]
        print(f'no {target} in {tree["data"]}')
        return None

    def check_exists(self):
        pass


if __name__ == "__main__":
    directory = '/home/mh491/Database'
    reader = BMRB_Reader(directory)
    reader.run()
