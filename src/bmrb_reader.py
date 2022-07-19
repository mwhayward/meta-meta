from pathlib import Path
import xml.etree.ElementTree as et
import nmrstarlib
import nmrpystar
import os
import sqlite3
from sqlite3 import Error


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
        metabolites = {}
        targets = {'name': ['save_assembly_1'],
                   'chemical_formula': ['Chem_comp.Formula'],
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
                for readfile in nmrstarlib.read_files(str(target_dir.joinpath(file))):
                    if 'NMR quality control of fragment libraries for screening' in readfile['save_entry_information'][
                       'Entry.Title']:
                        name = readfile['save_assembly_1']['Assembly.Name'].strip('\n')
                    else:
                        name = readfile['save_entry_information']['Entry.Title'].strip('\n')
                    if name not in metabolites.keys():
                        metabolites[name] = [readfile['data']]
                    else:
                        metabolites[name].append(readfile['data'])
        print(metabolites)
        print(count)

    def find_save_data(self, tree, target):
        # returns the first instance of save data that matches the given target
        for save in tree.keys():
            if not isinstance(tree[save], str):
                for entity in tree[save].keys():
                    if target == entity:
                        return tree[save][entity]
        print(f'no {target} in {tree["data"]}')

    def check_exists(self):
        pass


if __name__ == "__main__":
    directory = '/home/mh491/Database'
    reader = BMRB_Reader(directory)
    reader.run()
