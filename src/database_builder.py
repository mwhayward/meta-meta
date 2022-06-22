import pathlib
import statistics
import xml.etree.ElementTree as et
import sys
import os
import pandas as pd
import sqlite3
from sqlite3 import Error
import re


class Builder:
    def __init__(self, database):
        self.database = database
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
        self.normal_concentrations = None
        self.abnormal_concentrations = None

    def create_connection(self):
        # Creates a single connection to the database
        # This connection is used at the end of the script to export the tables to a db file
        conn = None
        try:
            conn = sqlite3.connect(self.database)
        except Error as e:
            print(e)
        return conn

    def build(self, count, elem):
        # A method that creates an entry for a single HMDB accession (from an etree element)
        # Calls the other methods in the Builder class to populate each table
        # Only gathers data from metabolites that match accession numbers to the files we have
        if elem.find('{http://www.hmdb.ca}accession').text in self.names.index:

            # Gather the bottom level metadata for the metabolite
            titles = self.metabolite_titles(elem)
            data = self.metabolite_data(count, elem)
            data = pd.DataFrame([data], columns=titles, index=[count])
            self.insert_into_table(data)

            # Gather the synonyms for this metabolite
            synonyms = self.get_synonyms(elem)
            for syn in synonyms:
                data = pd.DataFrame([[count, syn]], columns=['metabolite_id', 'synonyms'])
                if self.synonyms is None:
                    self.synonyms = data
                else:
                    self.synonyms = self.synonyms.append(data)

            # Gather ontology data for this metabolite
            ont = elem.find('{http://www.hmdb.ca}ontology')
            tag = '{http://www.hmdb.ca}term'
            self.retrieve_all(count, ont, tag)

            # Gather concentration data for this metabolite
            if elem.find('{http://www.hmdb.ca}normal_concentrations') is not None:
                concelem = elem.find('{http://www.hmdb.ca}normal_concentrations')
                for concentration in concelem.getchildren():
                    normconcs = self.get_concentrations(concentration, count)
                    if self.normal_concentrations is None:
                        self.normal_concentrations = normconcs
                    else:
                        self.normal_concentrations = self.normal_concentrations.append(normconcs)
            if elem.find('{http://www.hmdb.ca}abnormal_concentrations') is not None:
                concelem = elem.find('{http://www.hmdb.ca}abnormal_concentrations')
                for concentration in concelem.getchildren():
                    abnormconcs = self.get_concentrations(concentration, count)
                    if self.abnormal_concentrations is None:
                        self.abnormal_concentrations = abnormconcs
                    else:
                        self.abnormal_concentrations = self.abnormal_concentrations.append(abnormconcs)

            # Updates the id number and prints a report to the console
            print(f'parsed {count} metabolites out of 1430')
            count += 1
        return count

    def metabolite_titles(self, elem, parent=None):
        # A method for gathering column names for a table from a given etree element
        # Used for the bottom level metabolite metadata
        name = elem.tag.split('}', 1)[1]
        titles = ['metabolite_id']
        #titles = [f'{name}_id']
        for child in elem.getchildren():
            if len(child.getchildren()) == 0 and child.tag not in self.exceptions:
                name = child.tag.split('}', 1)[1]
                titles.append(name)
        if parent is not None:
            titles.append(f'{parent}_id')
        return titles

    def metabolite_data(self, id, elem):
        # A method for gathering values for a table from a given etree element
        # Used for the bottom level metabolite metadata
        data = [id]
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
    
    def get_concentrations(self, concelem, count):
        # A method for gathering concentration data
        # Used for both normal and abnormal concentration data
        # Does not add values to the concentration dataframe but returns a dataframe
        titles = ['metabolite_id']
        data = [count]
        for child in concelem.getchildren():
            if len(child.getchildren()) == 0:
                titles.append(child.tag.split('}', 1)[1])
                data.append(child.text)
        concs = pd.DataFrame([data], columns=titles)
        return concs

    def save_to_db(self):
        # Saves all dataframes to the db file
        # This method is called once from outside the class
        self.metabolites.to_sql('metabolites', self.conn, if_exists='replace', index=False)
        self.synonyms.to_sql('synonyms', self.conn, if_exists='replace', index=False)
        self.isin.to_sql('isin', self.conn, if_exists='replace', index=False)
        self.ontology.to_sql('ontology', self.conn, if_exists='replace', index=False)
        self.normal_concentrations.to_sql('normal_concentrations', self.conn, if_exists='replace', index=False)
        self.abnormal_concentrations.to_sql('abnormal_concentrations', self.conn, if_exists='replace', index=False)
        

def parse_metabolites(builder):
    # A method that parses through the single xml file and builds element trees for each metabolite
    # Calls the builder 'build' method to create a data entry for each metabolite
    # Contains a fixed count of 1430 metabolites as any parsing beyond this value results in sigkill
    directory = '/home/mh491/Database/HMDB_files'
    file = pathlib.Path(directory, 'hmdb_metabolites.xml')
    nsmap = {}
    count = 1
    for i, (event, elem) in enumerate(et.iterparse(file, events=('end', 'start-ns'))):
        if event == 'start-ns':
            ns, url = elem
            nsmap[ns] = url
        elif event == 'end' and elem.tag.split('}', 1)[1] in ['metabolite']:
            count = builder.build(count, elem)
            if count >= 1430:
                break

class Reader:
    def __init__(self, database):
        self.database = database
        self.conn = self.create_connection()
        self.spectra = None
        self.multiplets = None
        self.peaks = None

    def create_connection(self):
        # Creates a single connection to the database
        # This connection is used at the beginning of the script to collect accession numbers
        # and at the end of the script to export the tables to a db file
        conn = None
        try:
            conn = sqlite3.connect(self.database)
        except Error as e:
            print(e)
        return conn

    def run(self):
        # this method iterates through the accession numbers of meta data and attempts to gather chemical shift data
        # nmrML files are first priority, followed by text files then xml
        sql = 'select "metabolite_id", "accession" from metabolites'
        metabolites = pd.read_sql(sql, self.conn)
        nmrmlfiles = os.listdir('/home/mh491/Database/HMDB_files/nmrML_experimental_Feb15_2022')
        textfiles = os.listdir('/home/mh491/Database/HMDB_files/hmdb_nmr_peak_lists')
        xmlfiles = os.listdir('/home/mh491/Database/HMDB_files/xml_files')
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
                self.parsexml(xml_metabolites, metabolite_id)
    
    def parsenmrml(self, files, metabolite_id):
        directory = '/home/mh491/Database/HMDB_files/nmrML_experimental_Feb15_2022'
        for file in files:
            file = pathlib.Path(directory, file)
            if not file.parts[-1].endswith('.nmrML') or not '1H' in file.parts[-1]:
                continue
            try:
                tree = et.parse(file)
                root = tree.getroot()
                root_tag = root.tag.rstrip('nmrML')
                try:
                    frequency = next(root.iter(f'{root_tag}effectiveExcitationField')).get('value')
                except:
                    frequency = None
                try:
                    reference = next(root.iter(f'{root_tag}chemicalShiftStandard')).get('name')
                except:
                    reference = None
                ph = None
                titles = ['metabolite_id', 'frequency', 'reference', 'ph']
                spectrum_data = pd.DataFrame([[metabolite_id, frequency, reference, ph]], columns=titles)
                if self.spectra is None:
                    self.spectra = spectrum_data
                else:
                    self.spectra = self.spectra.append(spectrum_data)
                for i, multiplet in enumerate(root.iter(f'{root_tag}multiplet')):
                    multiplet_id = f'{metabolite_id}.{i + 1}'
                    center = multiplet.get('center')
                    atom_ref = multiplet.find(f'{root_tag}atoms').get('atomRefs')
                    multiplicity = multiplet.find(f'{root_tag}multiplicity').get('name')
                    titles = ['multiplet_id', 'metabolite_id', 'center', 'atom_ref', 'multiplicity']
                    multiplet_data = pd.DataFrame([[multiplet_id, metabolite_id, center, atom_ref, multiplicity]], columns=titles)
                    if self.multiplets is None:
                        self.multiplets = multiplet_data
                    else:
                        self.multiplets = self.multiplets.append(multiplet_data)
                    for j, peak in enumerate(multiplet.find(f'{root_tag}peakList').findall(f'{root_tag}peak')):
                        peak_id = f'{multiplet_id}.{j + 1}'
                        shift = peak.get('center')
                        intensity = peak.get('amplitude')
                        width = peak.get('width')
                        titles = ['peak_id', 'metabolite_id', 'multiplet_id', 'shift', 'intensity', 'width']
                        peak_data = pd.DataFrame([[peak_id, metabolite_id, multiplet_id, shift, intensity, width]], columns=titles)
                        if self.peaks is None:
                            self.peaks = peak_data
                        else:
                            self.peaks = self.peaks.append(peak_data)
            except Exception as e:
                print(f'Unable to parse file {file.name}')

    def parsetext(self, files, metabolite_id):
        # takes a set of text file filenames, gathers chemical shift data and formats it to fit the SQL schema
        directory = '/home/mh491/Database/HMDB_files/hmdb_nmr_peak_lists/'
        for file in files:
            file = directory+file
            locations = self.find_tables(file)
            multiplets = self.get_text_data(file, 'multiplets', locations)
            peaks = self.get_text_data(file, 'peaks', locations)
            if multiplets is not None:
                for index, multiplet in multiplets.iterrows():
                    ppmrange = [float(num) for num in multiplet['(ppm)'].split(' .. ')]
                    multiplet_id = f'{metabolite_id}.{index+1}'
                    if 'Shift1(ppm)' in multiplets.columns:
                        center = multiplet['Shift1(ppm)']
                    else:
                        center = statistics.mean(ppmrange)
                    if 'Atom1' in multiplets.columns:
                        atom_ref = multiplet['Atom1']
                    else:
                        atom_ref = None
                    multiplicity = multiplet['Type']
                    titles = ['multiplet_id', 'metabolite_id', 'center', 'atom_ref', 'multiplicity']
                    multiplet_data = pd.DataFrame([[multiplet_id, metabolite_id, center, atom_ref, multiplicity]], columns=titles)
                    if self.multiplets is None:
                        self.multiplets = multiplet_data
                    else:
                        self.multiplets = self.multiplets.append(multiplet_data)
                    for index, peak in peaks.iterrows():
                        if min(ppmrange) < float(peak[1]) < max(ppmrange):
                            peak_id = f'{multiplet_id}.{index + 1}'
                            shift = peak['(ppm)']
                            intensity = peak['Height']
                            width = None
                            titles = ['peak_id', 'metabolite_id', 'multiplet_id', 'shift', 'intensity', 'width']
                            peak_data = pd.DataFrame([[peak_id, metabolite_id, multiplet_id, shift, intensity, width]], columns=titles)
                            if self.peaks is None:
                                self.peaks = peak_data
                            else:
                                self.peaks = self.peaks.append(peak_data)

    def parsexml(self, files, metabolite_id):
        # method for parsing xml files from hmdb
        # by default, xml files do not have multiplet data so this method skips multiplets
        directory = '/home/mh491/Database/HMDB_files/xml_files'
        for file in files:
            file = pathlib.Path(directory, file)
            tree = et.parse(file)
            root = tree.getroot()
            frequency = self.get_element(root, 'frequency')[0]
            reference = self.get_element(root, 'chemical-shift-reference')[0]
            ph = self.get_element(root, 'sample-ph')[0]
            concentration = self.get_element(root, 'sample-concentration')[0]
            concentration_units = self.get_element(root, 'sample-concentration-units')[0]
            temperature = self.get_element(root, 'sample-temperature')[0]
            temperature_units = self.get_element(root, 'sample-temperature-units')[0]
            titles = ['metabolite_id', 'frequency', 'reference', 'ph', 'concentration', 'concentration_units', 'temperature', 'temperature_units']
            spectrum_data = pd.DataFrame([[metabolite_id, frequency, reference, ph, concentration, concentration_units, temperature, temperature_units]], columns=titles)
            if self.spectra is None:
                self.spectra = spectrum_data
            else:
                self.spectra = self.spectra.append(spectrum_data)
            for i, peak in enumerate(root.iter('nmr-one-d-peak')):
                peak_id = f'{metabolite_id}.1.{i}'
                shift = peak.find('chemical-shift').text
                intensity = peak.find('intensity').text
                width = None
                multiplet_id = None
                titles = ['peak_id', 'metabolite_id', 'multiplet_id', 'shift', 'intensity', 'width']
                peak_data = pd.DataFrame([[peak_id, metabolite_id, multiplet_id, shift, intensity, width]], columns=titles)
                if self.peaks is None:
                    self.peaks = peak_data
                else:
                    self.peaks = self.peaks.append(peak_data)

    def get_element(self, root, tag):
        out = []
        for elem in root.iter(tag):
            out.append(elem.text)
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


if __name__ == "__main__":
    directory = '/home/mh491/Database'
    database = f'{directory}/ccpn_metabolites.db'
    # builder = Builder(database)
    # parse_metabolites(builder)
    # builder.save_to_db()

    reader = Reader(database)
    reader.run()
    print(reader.spectra)
    print(reader.multiplets)
    print(reader.peaks)
    # nmrml_dir = directory+'/nmrML_experimental_Feb15_2022'
    # parse_hmdb_nmrml_data(directory, builder)
