from pathlib import Path
import statistics
import xml.etree.ElementTree as et
import sys
import os
import pandas as pd
import sqlite3
from sqlite3 import Error
import re


class Builder:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.conn = self.create_connection()
        self.exceptions = ['{http://www.hmdb.ca}secondary_accessions',
                           '{http://www.hmdb.ca}synonyms',
                           '{http://www.hmdb.ca}taxonomy',
                           '{http://www.hmdb.ca}ontology',
                           '{http://www.hmdb.ca}experimental_properties',
                           '{http://www.hmdb.ca}predicted_properties',
                           '{http://www.hmdb.ca}spectra',
                           '{http://www.hmdb.ca}biological_properties',
                           '{http://www.hmdb.ca}normal_concentrations',
                           '{http://www.hmdb.ca}abnormal_concentrations',
                           '{http://www.hmdb.ca}diseases',
                           '{http://www.hmdb.ca}general_references',
                           '{http://www.hmdb.ca}protein_associations']
        self.names = pd.read_csv('/home/mh491/Database/HMDB_files/metabolite_names.csv', index_col=0)
        self.metabolites = None
        self.synonyms = None
        self.isin = None
        self.ontology = None
        self.concentrations = None

    def create_connection(self):
        # Creates a single connection to the database
        # This connection is used at the end of the script to export the tables to a db file
        conn = None
        try:
            conn = sqlite3.connect(str(self.directory.joinpath('ccpn_metabolites.db')))
        except Error as e:
            print(e)
        return conn

    def build(self, count, elem):
        # A method that creates an entry for a single HMDB accession (from an etree element)
        # Calls the other methods in the Builder class to populate each table
        # Only gathers data from metabolites that match accession numbers to the files we have
        if elem.find('{http://www.hmdb.ca}accession').text in self.names.index:

            # Gather the bottom level metadata for the metabolite
            metabolite_id = f'SU:{count}'
            titles = self.metabolite_titles(elem)
            data = self.metabolite_data(metabolite_id, elem)
            data = pd.DataFrame([data], columns=titles, index=[metabolite_id])
            self.insert_into_table(data)

            # Gather the synonyms for this metabolite
            synonyms = self.get_synonyms(elem)
            for syn in synonyms:
                data = pd.DataFrame([[metabolite_id, syn]], columns=['metabolite_id', 'synonyms'])
                if self.synonyms is None:
                    self.synonyms = data
                else:
                    self.synonyms = self.synonyms.append(data)

            # Gather ontology data for this metabolite
            ont = elem.find('{http://www.hmdb.ca}ontology')
            tag = '{http://www.hmdb.ca}term'
            self.retrieve_all(metabolite_id, ont, tag)

            # Gather concentration data for this metabolite
            if elem.find('{http://www.hmdb.ca}normal_concentrations') is not None:
                concelem = elem.find('{http://www.hmdb.ca}normal_concentrations')
                for concentration in concelem.getchildren():
                    self.get_concentrations(concentration, metabolite_id)
            if elem.find('{http://www.hmdb.ca}abnormal_concentrations') is not None:
                concelem = elem.find('{http://www.hmdb.ca}abnormal_concentrations')
                for concentration in concelem.getchildren():
                    self.get_concentrations(concentration, metabolite_id)

            # Updates the id number and prints a report to the console
            print(f'parsed {count} metabolites out of 1430')
            count += 1
        return count

    def metabolite_titles(self, elem, parent=None):
        # A method for gathering column names for a table from a given etree element
        # Used for the bottom level metabolite metadata
        name = elem.tag.split('}', 1)[1]
        titles = ['metabolite_id']
        # titles = [f'{name}_id']
        for child in elem.getchildren():
            if len(child.getchildren()) == 0 and child.tag not in self.exceptions:
                name = child.tag.split('}', 1)[1]
                titles.append(name)
        if parent is not None:
            titles.append(f'{parent}_id')
        return titles

    def metabolite_data(self, metabolite_id, elem):
        # A method for gathering values for a table from a given etree element
        # Used for the bottom level metabolite metadata
        data = [metabolite_id]
        for child in elem.getchildren():
            if len(child.getchildren()) == 0 and child.tag not in self.exceptions:
                data.append(child.text)
        return data

    def insert_into_table(self, data):
        # A method for creating or appending data to the metabolite dataframe
        if self.metabolites is None:
            self.metabolites = data
        else:
            self.metabolites = self.metabolites.append(data)

    def get_synonyms(self, elem):
        # A method for gathering synonym data from a synonym element
        syn = elem.find('{http://www.hmdb.ca}synonyms')
        synonyms = []
        if len(syn.getchildren()) > 0:
            for child in syn.getchildren():
                synonyms.append(child.text)
        return synonyms

    def retrieve_all(self, metab_id, elem, tag=None, parent=None):
        # A method for gathering ontology data
        # Recursively scans the ontology element for nodes under the 'term' tag
        # Adds values to the ontology dataframe if they are not already present
        # Adds values to the isin dataframe
        if len(elem.getchildren()) > 0:
            for child in elem.getchildren():
                self.retrieve_all(metab_id, child, tag, elem)
        elif elem.tag == tag:
            titles = ['group', 'definition']
            data = pd.DataFrame([[elem.text, parent.find('{http://www.hmdb.ca}definition').text]], columns=titles)
            linktitles = ['metabolite_id', 'group']
            link = pd.DataFrame([[metab_id, elem.text]], columns=linktitles)
            if self.isin is None:
                self.isin = link
            else:
                self.isin = self.isin.append(link)
            if self.ontology is None:
                self.ontology = data
            elif elem.text not in self.ontology.index:
                self.ontology = self.ontology.append(data)

    def get_concentrations(self, concelem, metabolite_id):
        # A method for gathering concentration data
        # Used for both normal and abnormal concentration data
        # Does not add values to the concentration dataframe but returns a dataframe
        titles = ['metabolite_id']
        data = [metabolite_id]
        for child in concelem.getchildren():
            if len(child.getchildren()) == 0:
                titles.append(child.tag.split('}', 1)[1])
                data.append(child.text)
        concs = pd.DataFrame([data], columns=titles)
        if self.concentrations is None:
            self.concentrations = concs
        else:
            self.concentrations = self.concentrations.append(concs)

    def save_to_db(self):
        # Saves all dataframes to the db file
        # This method is called once from outside the class
        self.metabolites.to_sql('metabolites', self.conn, if_exists='replace', index=False)
        self.synonyms.to_sql('synonyms', self.conn, if_exists='replace', index=False)
        self.isin.to_sql('isin', self.conn, if_exists='replace', index=False)
        self.ontology.to_sql('ontology', self.conn, if_exists='replace', index=False)
        self.concentrations.to_sql('normal_concentrations', self.conn, if_exists='replace', index=False)

    def parse_metabolites(self):
        # A method that parses through the single xml file and builds element trees for each metabolite
        # Calls the builder 'build' method to create a data entry for each metabolite
        # Contains a fixed count of 1430 metabolites as any parsing beyond this value results in sigkill
        file = self.directory.joinpath('HMDB_files/hmdb_metabolites.xml')
        nsmap = {}
        count = 1
        for i, (event, elem) in enumerate(et.iterparse(file, events=('end', 'start-ns'))):
            if event == 'start-ns':
                ns, url = elem
                nsmap[ns] = url
            elif event == 'end' and elem.tag.split('}', 1)[1] in ['metabolite']:
                count = self.build(count, elem)
                if count >= 1430:
                    break


class Reader:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.conn = self.create_connection()

        self.sampletitles = ['sample_id', 'metabolite_id', 'pH', 'amount', 'amount_unit', 'reference']
        self.samples = pd.DataFrame(columns=self.sampletitles)
        self.sample_key = 1

        self.spectratitles = ['spectrum_id', 'sample_id', 'frequency', 'temperature', 'temperature_units']
        self.spectra = pd.DataFrame(columns=self.spectratitles)
        self.spectrum_key = 1

        self.multiplettitles = ['multiplet_id', 'spectrum_id', 'center', 'atom_ref', 'multiplicity']
        self.multiplets = pd.DataFrame(columns=self.multiplettitles)
        self.multiplet_key = 1

        self.peaktitles = ['peak_id', 'spectrum_id', 'multiplet_id', 'shift', 'intensity', 'width']
        self.peaks = pd.DataFrame(columns=self.peaktitles)
        self.peak_key = 1

    def create_connection(self):
        # Creates a single connection to the database
        # This connection is used at the beginning of the script to collect accession numbers
        # and at the end of the script to export the tables to a db file
        conn = None
        try:
            conn = sqlite3.connect(str(self.directory.joinpath('ccpn_metabolites.db')))
        except Error as e:
            print(e)
        return conn

    def run(self):
        # this method iterates through the accession numbers of meta data and attempts to gather chemical shift data
        # nmrML files are first priority, followed by text files then xml
        sql = 'select "metabolite_id", "accession" from metabolites'
        metabolites = pd.read_sql(sql, self.conn)
        nmrmlfiles = os.listdir(self.directory.joinpath('HMDB_files/nmrML_experimental_Feb15_2022'))
        textfiles = os.listdir(self.directory.joinpath('HMDB_files/hmdb_nmr_peak_lists'))
        xmlfiles = os.listdir(self.directory.joinpath('HMDB_files/xml_files'))
        for index, metabolite in metabolites.iterrows():
            metabolite_id = metabolite['metabolite_id']
            accession = metabolite['accession']
            nmrml_expr = re.compile(f'{accession}_.*_1H.nmrML')
            text_expr = re.compile(f'{accession}_nmroned.*')
            xml_expr = re.compile(f'{accession}_nmr_one_d.*')
            nmrml_metabolites = list(filter(nmrml_expr.match, nmrmlfiles))
            text_metabolites = list(filter(text_expr.match, textfiles))
            xml_metabolites = list(filter(xml_expr.match, xmlfiles))
            if len(nmrml_metabolites) > 0:
                self.parsenmrml(nmrml_metabolites, metabolite_id)
            elif len(text_metabolites) > 0:
                try:
                    self.parsetext(text_metabolites, metabolite_id)
                except:
                    print(f'error with {text_metabolites}')
            elif len(xml_metabolites) > 0:
                try:
                    self.parsexml(xml_metabolites, metabolite_id)
                except:
                    print(f'error with {xml_metabolites}')
            self.save_to_db()

    def parsenmrml(self, files, metabolite_id):
        # method for gathering sample/spectrum/multiplet/peak data from nmrML files
        # calls xml files to cover the shortcomings of nmrML with sample/spectrum data
        directory = self.directory.joinpath('HMDB_files/nmrML_experimental_Feb15_2022')
        for i, file in enumerate(files):
            sample_id = f'SA:{self.sample_key}'
            sample_data = {'sample_id': sample_id,
                           'metabolite_id': metabolite_id}
            spectrum_id = f'SP:{self.spectrum_key}'
            spectrum_data = {'spectrum_id': spectrum_id,
                             'sample_id': sample_id}
            file = directory.joinpath(file)
            if not str(file).endswith('.nmrML') or not '1H' in str(file):
                continue
            tree = et.parse(file)
            root = tree.getroot()
            root_tag = root.tag.rstrip('nmrML')
            try:
                sample_data['reference'] = next(root.iter(f'{root_tag}chemicalShiftStandard')).get('name')
            except:
                pass
            try:
                spectrum_data['frequency'] = next(root.iter(f'{root_tag}effectiveExcitationField')).get('value')
            except:
                pass
            self.get_xml_spectrum_data(file, sample_data, spectrum_data)
            peak_count = 1
            for j, multiplet in enumerate(root.iter(f'{root_tag}multiplet')):
                multiplet_id = f'MT:{spectrum_id.split(":")[-1]}.{j + 1}'
                multiplet_data = {'multiplet_id': multiplet_id,
                                  'spectrum_id': spectrum_id,
                                  'center': multiplet.get('center'),
                                  'atom_ref': multiplet.find(f'{root_tag}atoms').get('atomRefs'),
                                  'multiplicity': multiplet.find(f'{root_tag}multiplicity').get('name')}
                self.multiplets.loc[self.multiplet_key] = multiplet_data
                self.multiplet_key += 1
                for k, peak in enumerate(multiplet.find(f'{root_tag}peakList').findall(f'{root_tag}peak')):
                    peak_data = {'peak_id': f'PK{spectrum_id.split(":")[-1]}.{peak_count + 1}',
                                 'spectrum_id': spectrum_id,
                                 'multiplet_id': multiplet_id,
                                 'shift': peak.get('center'),
                                 'intensity': peak.get('amplitude'),
                                 'width': peak.get('width')}
                    self.peaks.loc[self.peak_key] = peak_data
                    self.peak_key += 1
                    peak_count += 1

    def parsetext(self, files, metabolite_id):
        # takes a set of text file filenames, gathers chemical shift data and formats it to fit the SQL schema
        directory = self.directory.joinpath('HMDB_files/hmdb_nmr_peak_lists')
        for i, file in enumerate(files):
            sample_id = f'SA:{self.sample_key}'
            spectrum_id = f'SP:{self.spectrum_key}'
            sample_data = {'sample_id': sample_id}
            spectrum_data = {'spectrum_id': spectrum_id,
                             'sample_id': sample_id}
            file = directory.joinpath(file)
            self.get_xml_spectrum_data(file, sample_data, spectrum_data)
            locations = self.find_tables(file)
            multiplets = self.get_text_data(file, 'multiplets', locations)
            peaks = self.get_text_data(file, 'peaks', locations)
            peak_count = 1
            if multiplets is not None:
                for j, multiplet in multiplets.iterrows():
                    ppmrange = [float(num) for num in multiplet['(ppm)'].split(' .. ')]
                    multiplet_id = f'MT:{spectrum_id.split(":")[-1]}.{j + 1}'
                    multiplet_data = {'multiplet_id': multiplet_id,
                                      'spectrum_id': spectrum_id,}
                    if 'Shift1(ppm)' in multiplets.columns:
                        multiplet_data['center'] = multiplet['Shift1(ppm)']
                    else:
                        multiplet_data['center'] = statistics.mean(ppmrange)
                    if 'Atom1' in multiplets.columns:
                        multiplet_data['atom_ref'] = multiplet['Atom1']
                    multiplet_data['multiplicity'] = multiplet['Type']
                    self.multiplets.loc[self.multiplet_key] = multiplet_data
                    self.multiplet_key += 1
                    for k, peak in peaks.iterrows():
                        if min(ppmrange) < float(peak[1]) < max(ppmrange):
                            peak_data = {'peak_id': f'PK{spectrum_id.split(":")[-1]}.{peak_count + 1}',
                                         'spectrum_id': spectrum_id,
                                         'multiplet_id': multiplet_id,
                                         'shift': peak['(ppm)'],
                                         'intensity': peak['Height']}
                            self.peaks.loc[self.peak_key] = peak_data
                            self.peak_key += 1
                            peak_count += 1

    def parsexml(self, files, metabolite_id):
        # method for parsing xml files from hmdb
        # by default, xml files do not have multiplet data so this method skips multiplets
        for i, file in enumerate(files):
            file = self.directory.joinpath(f'HMDB_files/xml_files/{file}')
            sample_id = f'SA:{self.sample_key}'
            spectrum_id = f'SP:{self.spectrum_key}'
            sample_data = {'sample_id': sample_id}
            spectrum_data = {'spectrum_id': spectrum_id,
                             'sample_id': sample_id}
            root = self.get_xml_spectrum_data(file, sample_data, spectrum_data)
            for j, peak in enumerate(root.iter('nmr-one-d-peak')):
                peak_data = {'peak_id': f'PK{spectrum_id.split(":")[-1]}.{j + 1}',
                             'spectrum_id': spectrum_id,
                             'shift': peak.find('chemical-shift').text,
                             'intensity': peak.find('intensity').text}
                self.peaks.loc[self.peak_key] = peak_data
                self.peak_key += 1

    def get_xml_spectrum_data(self, file, sample_data, spectrum_data):
        # in a separate method to allow all file reader methods to check metadata from the relevant xml file
        # tries to get xml experiment data but prioritises already existing data from nmrML files
        if file.suffix == '.xml':
            specnum = None
        if file.suffix == '.nmrML':
            specnum = str(file).split('_')[-2]
        elif file.suffix == '.txt':
            specnum = str(file).split('_')[-2]
        if specnum is not None:
            xmlfiles = os.listdir(self.directory.joinpath('HMDB_files/xml_files'))
            filetargets = [f for f in xmlfiles if f'one_d_spectrum_{specnum}' in f]
            if len(filetargets) > 0:
                file = self.directory.joinpath(f'HMDB_files/xml_files/{filetargets[0]}')
            else:
                # print(f'no xml file to complement {file}')
                self.samples.loc[self.sample_key] = sample_data
                self.sample_key += 1
                self.spectra.loc[self.spectrum_key] = spectrum_data
                self.spectrum_key += 1
                return None
        tree = et.parse(file)
        root = tree.getroot()
        xml_sample = {'pH': self.get_element(root, 'sample-ph')[0],
                      'amount': self.get_element(root, 'sample-concentration')[0],
                      'amount_units': self.get_element(root, 'sample-concentration-units')[0],
                      'reference': self.get_element(root, 'chemical-shift-reference')[0]}
        xml_spectrum = {'temperature': self.get_element(root, 'sample-temperature')[0],
                        'temperature_units': self.get_element(root, 'sample-temperature-units')[0],
                        'frequency': self.get_element(root, 'frequency')[0]}
        sample_data.update(xml_sample)
        spectrum_data.update(xml_spectrum)
        self.samples.loc[self.sample_key] = sample_data
        self.sample_key += 1
        self.spectra.loc[self.spectrum_key] = spectrum_data
        self.spectrum_key += 1
        return root

    def get_element(self, root, tag):
        # retrieves data from an element in the et tree
        # or returns None if it doesn't exist
        out = []
        for elem in root.iter(tag):
            out.append(elem.text)
        if len(out) == 0:
            out = [None]
        return out

    def get_text_data(self, file, feature, locations):
        # gathers either peak or multiplet data from text files
        # returns a pandas dataframe identical to the one from the text file
        if locations is 'C13':
            return None
        startline = locations[feature]
        if startline == -1:
            print(f'no {feature} table in {file}')
            return None
        file = open(file, 'r')
        titles = []
        table = []
        start = False
        for i, line in enumerate(file.readlines()):
            line = line.lstrip(' ')
            if i == startline:
                start = True
            if line != '\n':
                line = line.rstrip('\n')
            if line.startswith('No') and start is True:
                for value in line.split('\t'):
                    titles.append(value.replace(' ', ''))
                while '' in titles:
                    titles.remove('')
            elif line[0].isdigit() and start is True:
                row = line.split('\t')
                while '' in row:
                    row.remove('')
                table.append(row)
            elif start is True and bool(table) is not False:
                break
        df = pd.DataFrame(table, columns=titles)
        if df.empty:
            print(f'empty multiplet table in file: {file}')
            return None
        return df

    def find_tables(self, file):
        # finds the line where the peak and multiplet tables start in the text file
        # returns a dictionary on success or a string on fail
        locations = {}
        file = open(file, 'r')
        for i, line in enumerate(file.readlines()):
            if line.startswith('DUoptxwinnmr'):
                return 'C13'
            if 'peaks' in line.casefold():
                locations['peaks'] = i
            if 'multiplets' in line.casefold():
                locations['multiplets'] = i
        if 'peaks' not in locations:
            locations['peaks'] = 0
        if 'multiplets' not in locations:
            locations['multiplets'] = -1
        return locations

    def save_to_db(self):
        # saves the tables generated by the reader to the db file
        self.samples.to_sql('samples', self.conn, if_exists='replace', index=False)
        self.spectra.to_sql('spectra', self.conn, if_exists='replace', index=False)
        self.multiplets.to_sql('multiplets', self.conn, if_exists='replace', index=False)
        self.peaks.to_sql('peaks', self.conn, if_exists='replace', index=False)


if __name__ == "__main__":
    directory = '/home/mh491/Database'
    # builder = Builder(directory)
    # builder.parse_metabolites()
    # builder.save_to_db()

    reader = Reader(directory)
    reader.run()
    print(reader.samples)
    print(reader.spectra)
    print(reader.multiplets)
    print(reader.peaks)
    # nmrml_dir = directory+'/nmrML_experimental_Feb15_2022'
    # parse_hmdb_nmrml_data(directory, builder)
