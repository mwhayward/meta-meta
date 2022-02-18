import pathlib
import xml.etree.ElementTree as et
import sqlite3
from sqlite3 import Error

class HMDB_nmrML_Reader:

    def __init__(self, directory):
        self._directory = directory
        self.conn = None

    def create_database(self):
        sql_create_metabolites_table = """ CREATE TABLE IF NOT EXISTS metabolites (
                                                        metabolite_id integer PRIMARY KEY,
                                                        hmdb_id text NOT NULL,
                                                        frequency float,
                                                        reference text
                                                    );"""

        sql_create_multiplets_table = """CREATE TABLE IF NOT EXISTS multiplets (
                                                    multiplet_id text PRIMARY KEY,
                                                    metabolite_id integer NOT NULL,
                                                    center float NOT NULL,
                                                    atom_ref integer,
                                                    multiplicity text, 
                                                    FOREIGN KEY (metabolite_id) REFERENCES metabolites (metabolite_id)
                                                );"""

        sql_create_peaks_table = """CREATE TABLE IF NOT EXISTS peaks (
                                                    peak_id text PRIMARY KEY,
                                                    metabolite_id integer NOT NULL,
                                                    multiplet_id text NOT NULL,
                                                    shift float NOT NULL,
                                                    intensity float NOT NULL,
                                                    width float,
                                                    FOREIGN KEY (metabolite_id) REFERENCES metabolites (metabolite_id)
                                                    FOREIGN KEY (multiplet_id) REFERENCES multiplets (multiplet_id)
                                                );"""
        if self.conn is not None:
            self.create_table(sql_create_metabolites_table)
            self.create_table(sql_create_multiplets_table)
            self.create_table(sql_create_peaks_table)
        else:
            print("Error! cannot create the database connection.")

    def create_table(self, create_table_sql):
        try:
            c = self.conn.cursor()
            c.execute(create_table_sql)
        except Error as e:
            print(e)

    def metabolite_output(self, metabolite_data):
        sql = ''' INSERT INTO metabolites(metabolite_id, hmdb_id, frequency, reference)
                  VALUES(?,?,?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, metabolite_data)
        self.conn.commit()
        return cur.lastrowid

    def multiplet_output(self, multiplet_data):
        sql = ''' INSERT INTO multiplets(multiplet_id, metabolite_id, center, atom_ref, multiplicity)
                          VALUES(?,?,?,?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, multiplet_data)
        self.conn.commit()
        return cur.lastrowid

    def peak_output(self, peak_data):
        sql = ''' INSERT INTO peaks(peak_id, metabolite_id, multiplet_id, shift, intensity, width)
                          VALUES(?,?,?,?,?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, peak_data)
        self.conn.commit()
        return cur.lastrowid

    def run(self):
        directory = pathlib.Path(self._directory)
        database = f'{directory}/hmdb_nmrML_metabolites.db'
        print(database)
        try:
            self.conn = sqlite3.connect(database)
        except Error as e:
            print(e)
        self.create_database()
        
        for file in directory.iterdir():
            # todo: reintroduce screener statement below
            if not file.parts[-1].endswith('.nmrML') or not '1H' in file.parts[-1]:
                continue

            try:
                tree = et.parse(file)
            except Exception as e:
                print(f'Unable to parse file {file.name}')
                pass

            metabolite_id = file.stem.split('_')[1]
            hmdb_id = file.stem.split('_')[0]
            root = tree.getroot()
            frequency = None
            reference = None
            feild_branch = '{http://nmrml.org/schema}acquisition/' \
                           '{http://nmrml.org/schema}acquisition1D/' \
                           '{http://nmrml.org/schema}acquisitionParameterSet/' \
                           '{http://nmrml.org/schema}DirectDimensionParameterSet/' \
                           '{http://nmrml.org/schema}effectiveExcitationField'
            if root.find(feild_branch) != None:
                frequency = (root.find(feild_branch).get('value'))
            reference_branch = '{http://nmrml.org/schema}sampleList/ \
                                {http://nmrml.org/schema}sample/ \
                                {http://nmrml.org/schema}chemicalShiftStandard'
            if root.find(reference_branch) != None:
                reference = (root.find(reference_branch).get('name'))
            metabolite_data = (metabolite_id, hmdb_id, frequency, reference)
            self.metabolite_output(metabolite_data)
            for i, multiplet in enumerate(root.iter('multiplet')):
                multiplet_id = f'{metabolite_id}.{i+1}'
                center = multiplet.get('center')
                atom_ref = multiplet.find('atoms').get('atomRefs')
                multiplicity = multiplet.find('multiplicity').get('name')
                multiplet_data =(multiplet_id, metabolite_id, center, atom_ref, multiplicity)
                try:
                    self.multiplet_output(multiplet_data)
                except:
                    print(f'error with multiplet {multiplet_id} in file {file.name}')
                for j, peak in enumerate(multiplet.find('peakList').findall('peak')):
                    peak_id = f'{multiplet_id}.{j+1}'
                    shift = peak.get('center')
                    intensity = peak.get('amplitude')
                    width = peak.get('width')
                    peak_data = (peak_id, metabolite_id, multiplet_id, shift, intensity, width)
                    try:
                        self.peak_output(peak_data)
                    except:
                        print(f'error with peak {peak_id} in file {file.name}')
                    



if __name__ == "__main__":
    # directory = '/home/mh491/Documents/nmrML_practice folder'
    directory = '/home/mh491/Documents/nmrML_experimental_Feb15_2022'
    reader = HMDB_nmrML_Reader(directory)
    reader.run()