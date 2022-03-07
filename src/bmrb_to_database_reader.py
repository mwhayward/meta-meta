import pathlib
import xml.etree.ElementTree as et
import nmrstarlib
import sys
import statistics
import sqlite3
from sqlite3 import Error


class BMRB_Reader:

    def __init__(self, directory):
        self._directory = directory
        self.conn = None

    def has_1h(self, tree):
        result = False
        SAVE_FRAME = 'save_entry_information'
        for top_frame in [top for top in tree if top.startswith(SAVE_FRAME)]:
            for loop in [elem for elem in tree[top_frame] if elem.startswith('loop')]:
                for elem in tree[top_frame][loop][1:]:
                    for i, frame in enumerate(elem):
                        for key in frame.keys():
                            if key == 'Datum.Type':
                                if frame[key] == '1H chemical shifts' and result != True:
                                    result = True
        return result

    def get_1h_shifts(self, tree):
        result = []
        SAVE_FRAME = 'save_spectral_peak_1H'
        for top_frame in [top for top in tree if top == SAVE_FRAME]:
            for loop in [elem for elem in tree[top_frame] if elem.startswith('loop')]:
                for elem in tree[top_frame][loop][1:]:
                    for i, line in enumerate(elem):
                        for key in line.keys():
                            if key == 'Spectral_transition_char.Chem_shift_val':
                                try:
                                    result.append(float(line[key]))
                                except:
                                    print(f"   couldn\'t convert {line['Atom_chem_shift.Val']} to float")
        return result

    def drop_tables(self):
        try:
            c = self.conn.cursor()
            c.execute('DROP TABLE metabolites')
            c.execute('DROP TABLE peaks')
        except Error as e:
            print(e)

    def create_database(self):
        sql_create_metabolites_table = """ CREATE TABLE IF NOT EXISTS metabolites (
                                                        bmrb_id integer PRIMARY KEY,
                                                        metabolite_name text
                                                    );"""

        sql_create_multiplets_table = """CREATE TABLE IF NOT EXISTS multiplets (
                                                    multiplet_id text PRIMARY KEY,
                                                    bmrb_id integer NOT NULL,
                                                    center float NOT NULL,
                                                    atom_ref integer,
                                                    multiplicity text, 
                                                    FOREIGN KEY (bmrb_id) REFERENCES metabolites (metabolite_id)
                                                );"""

        sql_create_peaks_table = """CREATE TABLE IF NOT EXISTS peaks (
                                                    peak_id text PRIMARY KEY,
                                                    bmrb_id integer NOT NULL,
                                                    shift float NOT NULL,
                                                    FOREIGN KEY (bmrb_id) REFERENCES multiplets (multiplet_id)
                                                );"""
        if self.conn is not None:
            self.create_table(sql_create_metabolites_table)
            #self.create_table(sql_create_multiplets_table)
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
        sql = ''' INSERT INTO metabolites(bmrb_id, metabolite_name)
                  VALUES(?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, metabolite_data)
        self.conn.commit()
        return cur.lastrowid

    def multiplet_output(self, multiplet_data):
        sql = ''' INSERT INTO multiplets(multiplet_id, bmrb_id, center, atom_ref, multiplicity)
                          VALUES(?,?,?,?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, multiplet_data)
        self.conn.commit()
        return cur.lastrowid

    def peak_output(self, peak_data):
        sql = ''' INSERT INTO peaks(peak_id, bmrb_id, shift)
                          VALUES(?,?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, peak_data)
        self.conn.commit()
        return cur.lastrowid

    def id_to_index(self, id):
        return id[4:].lstrip('0')

    def run(self):
        directory = pathlib.Path(self._directory)
        database = f'{directory}/bmrb_str_metabolites.db'
        try:
            self.conn = sqlite3.connect(database)
        except Error as e:
            print(e)
        self.drop_tables()
        self.create_database()
        print(f"processing files in {directory}")

        file_ids = [int(self.id_to_index(str(file.stem))) for file in directory.iterdir() if file.suffix == '.str']
        for id in sorted(file_ids):
            file = pathlib.Path(directory, f'bmse{str(id).zfill(6)}.str')
            tree = next(nmrstarlib.read_files(str(file)))
            name = tree['save_entry_information']['Entry.Title'].replace('_', '-').strip('\n')
            shifts = self.get_1h_shifts(tree)
            molecule_id = self.id_to_index(tree.id)

            if len(shifts) == 0:
                print(f"ignoring file {file} it doesn't contain any assigned 1h shifts")
                continue

            if self.has_1h(tree):
                metabolite_data = (molecule_id, name)
                self.metabolite_output(metabolite_data)
                for num, peak in enumerate(shifts):
                    peak_id = f'{molecule_id}.{num+1}'
                    peak_data = (peak_id, molecule_id, peak)
                    self.peak_output(peak_data)

if __name__ == "__main__":
    directory = '/home/mh491/Metameta_Files/bmrb_nmr_spectra'
    reader = BMRB_Reader(directory)
    reader.run()
