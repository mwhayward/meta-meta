import pathlib
import xml.etree.ElementTree as et
import nmrstarlib
import sys
import statistics
import sqlite3
from sqlite3 import Error

class HMDB_Metabolite_Reader:

    def __init__(self, directory):
        self._directory = directory

    def fix_tag(self, ns, tag, nsmap):
        return '{' + nsmap[ns] + '}' + tag

    def get_tag(self, elem):
        return elem.tag.split('}', 1)[1]

    def get_text(self, file, tag):
        result = []
        for i, (event, elem) in enumerate(et.iterparse(file, events=('end', 'tail', 'text'))):
            if event == 'text':
                result.append("".join(elem.itertext()))
            elif event == 'end':
                elem.clear()
                return ''.join(result)

    def id_to_index(self, id):
        return id[4:].lstrip('0')

    def output(self, id, name):
        print(f'HMDB-{self.id_to_index(id)},{name}')

    def run(self):
        file = pathlib.Path(self._directory, 'hmdb_metabolites.xml')

        nsmap = {}

        for i, (event, elem) in enumerate(et.iterparse(file, events=('end', 'start-ns'))):
            if event == 'start-ns':
                ns, url = elem
                nsmap[ns] = url
            elif event == 'end' and self.get_tag(elem) in ['metabolite']:
                name = None
                accession = None
                for accession in elem.findall('{http://www.hmdb.ca}accession'):
                    accession = accession.text
                for name in elem.findall('{http://www.hmdb.ca}name'):
                    name = name.text.replace(',', '-')
                self.output(accession, name)
                elem.clear()


class HMDB_Metabolites_to_CSV(HMDB_Metabolite_Reader):
    def __init__(self, directory):
        super(HMDB_Metabolites_to_CSV, self).__init__(directory)
        self.out_file = 'id_name.csv'

    def run(self):
        with open(pathlib.Path(directory, self.out_file), 'w') as self.csv_file:
            super(HMDB_Metabolites_to_CSV, self).run()

    def output(self, id, name):
        print(f'{id},{name}', file=self.csv_file)


class HMDB_Reader:

    def __init__(self, directory):
        self._directory = directory

    def shifts_as_list(self, peak_list):
        return ', '.join(["%7.3f" % peak[0] for peak in peak_list])

    def intensities_as_list(self, peak_list):
        return ', '.join(["%7.3f" % peak[-1] for peak in peak_list])

    def output(self, file, spectrum_id, molecule_id, frequency, ph, reference, temperature,multiplet_list, peak_list):
        print(f"file: {file}")
        print(f"id:   HMDB-{molecule_id}")
        shifts = self.shifts_as_list(peak_list)
        print(f"shifts: {shifts}")
        print()

    def run(self):
        directory = pathlib.Path(self._directory)
        print(f"processing files in {directory}")

        for file in directory.iterdir():
            if not file.parts[-1].endswith('.xml') or not 'nmr_one_d_spectrum' in file.parts[-1]:
                continue

            tree = et.parse(file)
            root = tree.getroot()
            for spectrum in root.iter('nmr-one-d'):
                # TODO: add variable to read nucleus of interest eg. 1H, 13C etc
                if spectrum.find('nucleus').text == "1H":
                    spectrum_id = spectrum.find('id').text
                    molecule_id = spectrum.find('database-id').text
                    if spectrum.find('frequency').text:
                        frequency = spectrum.find('frequency').text.rstrip(' MHz')
                    else:
                        frequency = None
                    ph = spectrum.find('sample-ph').text
                    reference = spectrum.find('chemical-shift-reference').text
                    temperature = spectrum.find('sample-temperature').text
                    peak_list = []
                    multiplet_list = []

                    multiplets = spectrum.findall("nmr-one-d-peaks")
                    for i, multiplet in enumerate(multiplets):
                        midpoint_list = []
                        for j, peak in enumerate(multiplet):
                            peak_id = f'{spectrum_id}.{i+1}.{j+1}'
                            shift = float(peak.find("chemical-shift").text)
                            if peak.find("intensity").text:
                                intensity = float(peak.find("intensity").text)
                            else:
                                intensity = None
                            peak_list.append((peak_id, shift, intensity))
                            midpoint_list.append(shift)
                        midpoint = statistics.median(midpoint_list)
                        multiplet_list.append((f'{spectrum_id}.{i+1}', midpoint))

                    # peak_list = sorted(peak_list)

                    self.output(file, spectrum_id, molecule_id, frequency, ph, reference, temperature, multiplet_list, peak_list)


class HMDB_to_MYSQL(HMDB_Reader):

    def __init__(self, directory):
        super(HMDB_to_MYSQL, self).__init__(directory)
        self.conn = None
        # self.metabolites_out_file = 'hmdb_xml_metbolites.csv'
        # self. multiplets_out_file = 'hmdb_xml_multiplets.csv'
        # self.peaks_out_file = 'hmdb_xml_peaks.csv'

    def create_database(self):
        sql_create_metabolites_table = """ CREATE TABLE IF NOT EXISTS metabolites (
                                                        metabolite_id integer PRIMARY KEY,
                                                        hmdb_id text NOT NULL,
                                                        field integer,
                                                        ph float,
                                                        reference text,
                                                        temperature integer
                                                    ); """

        sql_create_multiplets_table = """CREATE TABLE IF NOT EXISTS multiplets (
                                                    multiplet_id float PRIMARY KEY,
                                                    metabolite_id integer NOT NULL,
                                                    center float NOT NULL,
                                                    FOREIGN KEY (metabolite_id) REFERENCES metabolites (id)
                                                );"""

        sql_create_peaks_table = """CREATE TABLE IF NOT EXISTS peaks (
                                                    peak_id float PRIMARY KEY,
                                                    metabolite_id integer NOT NULL,
                                                    shift float NOT NULL,
                                                    intensity float NOT NULL,
                                                    FOREIGN KEY (metabolite_id) REFERENCES metabolites (id)
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

    def run(self):
        database = r'/home/mh491/Metameta_Files/hmdb_nmr_spectra/hmdb_metabolites.db'
        try:
            self.conn = sqlite3.connect(database)
        except Error as e:
            print(e)
        self.create_database()

        #with open(pathlib.Path(directory, self.metabolites_out_file), 'w') as self.metabolites_csv_file, \
        #     open(pathlib.Path(directory, self.multiplets_out_file), 'w') as self.multiplets_csv_file, \
        #     open(pathlib.Path(directory, self.peaks_out_file), 'w') as self.peaks_csv_file:
        super(HMDB_to_MYSQL, self).run()

    def output(self, file, spectrum_id, molecule_id, frequency, ph, reference, temperature, multiplet_list, peak_list):
        # shifts_text = self.shifts_as_list(peak_list)
        # intensities_text = self.intensities_as_list(peak_list)
        # print(f' {spectrum_id}, {molecule_id}, {frequency}, {ph}, {reference}, {temperature}', file=self.metabolites_csv_file)
        # for multiplet in multiplet_list:
        #     print(f'{multiplet[0]}, {spectrum_id}, {multiplet[-1]}', file=self.multiplets_csv_file)
        # for peak in peak_list:
        #     print(f'{peak[0]}, {spectrum_id}, {peak[1]}, {peak[-1]}', file=self.peaks_csv_file)
        metabolite_data = (spectrum_id, molecule_id, frequency, ph, reference, temperature)
        try:
            self.create_metabolite(metabolite_data)
        except Error as e:
            print(e)
            print(f'error with metabolite: {spectrum_id}')
        for multiplet in multiplet_list:
            multiplet_data = (multiplet[0], spectrum_id, multiplet[-1])
            try:
                self.create_multiplet(multiplet_data)
            except Error as e:
                print(e)
                print(f'error with multiplet: {multiplet_data[0]}')
        for peak in peak_list:
            peak_data = (peak[0], spectrum_id, peak[1], peak[-1])
            try:
                self.create_peak(peak_data)
            except Error as e:
                print(e)
                print(f'error with peak: {peak_data[0]}')

    def create_metabolite(self, metabolite):
        sql = ''' INSERT INTO metabolites(metabolite_id, hmdb_id, field, ph, reference, temperature)
                  VALUES(?,?,?,?,?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, metabolite)
        self.conn.commit()
        return cur.lastrowid

    def create_multiplet(self, multiplet):
        sql = ''' INSERT INTO multiplets(multiplet_id, metabolite_id, center)
                  VALUES(?,?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, multiplet)
        self.conn.commit()
        return cur.lastrowid

    def create_peak(self, peak):
        sql = ''' INSERT INTO peaks(peak_id, metabolite_id, shift, intensity)
                  VALUES(?,?,?,?) '''
        cur = self.conn.cursor()
        cur.execute(sql, peak)
        self.conn.commit()
        return cur.lastrowid


class BMRB_Reader:

    def __init__(self, directory):
        self._directory = directory

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
        for top_frame in [top for top in tree if top.startswith(SAVE_FRAME)]:
            for loop in [elem for elem in tree[top_frame] if elem.startswith('loop')]:
                for elem in tree[top_frame][loop][1:]:
                    for i, line in enumerate(elem):
                        for key in line.keys():
                            if key == 'Spectral_transition_char.Chem_shift_val':
                                try:
                                    result.append(float(line['Spectral_transition_char.Chem_shift_val']))
                                except:
                                    print(f"   couldn\'t convert {line['Atom_chem_shift.Val']} to float")
        return result

    def shifts_as_list(self, shifts):
        return ', '.join(["%7.3f" % shift for shift in sorted(shifts)])

    def output(self, file, moleule_id, spectrum_id, shifts, name):
        print(f"file: {file}")
        print(f"  id: BMRB-{id}")
        print(f"  shifts: {self.shifts_as_list(shifts)}")
        print(f"  molecule: {name}")
        print()

    def id_to_index(self, id):
        return id[4:].lstrip('0')

    def run(self):
        directory = pathlib.Path(self._directory)
        print(f"processing files in {directory}")

        file_ids = [int(self.id_to_index(str(file.stem))) for file in directory.iterdir() if file.suffix == '.str']
        for id in sorted(file_ids):
            file = pathlib.Path(directory, f'bmse{str(id).zfill(6)}.str')
            tree = next(nmrstarlib.read_files(str(file)))

            name = tree['save_entry_information']['Entry.Title']
            shifts = self.get_1h_shifts(tree)
            molecule_id = self.id_to_index(tree.id)

            if len(shifts) == 0:
                print(f"ignoring file {file} it doesn't contain any assigned 13c shifts")
                continue

            if self.has_1h(tree):
                self.output(file, molecule_id, 1, shifts, name)


# TODO: deal with multiple frames report peak shifts even if not assigned
class BMRB_to_CSV(BMRB_Reader):

    def __init__(self, directory):
        super(BMRB_to_CSV, self).__init__(directory)
        self.id_shifts_file = pathlib.Path(directory, 'id_shifts.csv')
        self.id_name_file = pathlib.Path(directory, 'id_name.csv')

    def run(self):
        with open(self.id_shifts_file, 'w') as self.id_shifts_csv_file, \
                open(self.id_name_file, 'w') as self.id_name_csv_file:
            super(BMRB_to_CSV, self).run()

    def output(self, file, molecule_id, spectrum_id, shifts, name):
        # shifts_text = self.shifts_as_list(shifts)
        name = name.replace('_', ' ')
        print(f'BMRB-{molecule_id}, {spectrum_id}, {self.shifts_as_list(shifts)}', file=self.id_shifts_csv_file)
        print(f'BMRB-{molecule_id}, "{name}"', file=self.id_name_csv_file)


class MMCD_Parser:
    def __init__(self, file):
        self._file = file
        self.shifts = None
        self.name = None

    def parse(self):
        result = []
        format = 2
        with open(self._file) as mmcd_file:
            for i, line in enumerate(mmcd_file):
                if i == 0:
                    for elem in line.strip().split(','):
                        if elem.strip().startswith('NAME'):
                            self.name = elem.split('=')[1].strip()
                        if elem.strip().startswith('DU'):
                            format = 1

                if format == 1 and i > 3:
                    fields = line.strip().split()
                    if len(fields) > 0:
                        result.append(float(fields[3]))
                elif format == 2 and i > 1:
                    fields = line.strip().split()
                    if len(fields) > 0:
                        result.append(float(fields[1]))

        self.shifts = sorted(result)
        return result


class MMCD_Reader:

    def __init__(self, directory):
        self._directory = directory

    def output(self, file, molecule_id, spectrum_id, shifts, name):
        print(f"file: {file}")
        print(f"  id: MMCD-{molecule_id}")
        print(f'  name: {name}')
        parser = MMCD_Parser(file)
        parser.parse()

        print(f"  shifts: {parser.shifts}")

    def run(self):
        directory = pathlib.Path(self._directory)
        print(f"processing files in {directory}")

        for i, file in enumerate(sorted(directory.iterdir())):
            if not file.stem.startswith('expnmr_'):
                continue

            molecule_id = file.stem.split('_')[1].lstrip('0')
            parser = MMCD_Parser(file)
            parser.parse()

            self.output(file, molecule_id, 1, parser.shifts, parser.name)


class MMCD_to_CSV(MMCD_Reader):

    def __init__(self, directory):
        super(MMCD_to_CSV, self).__init__(directory)
        self.id_shifts_file = pathlib.Path(directory, 'id_shifts.csv')
        self.id_name_file = pathlib.Path(directory, 'id_name.csv')

    def run(self):
        with open(self.id_shifts_file, 'w') as self.id_shifts_csv_file, \
                open(self.id_name_file, 'w') as self.id_name_csv_file:
            super(MMCD_to_CSV, self).run()

    def shifts_as_list(self, shifts):
        return ', '.join(["%7.3f" % shift for shift in shifts])

    def output(self, file,molecule_id, spectrum_id, shifts, name):
        # shifts_text = self.shifts_as_list(shifts)
        if not (name.startswith('expnmr_') or name.startswith('cq_')):
            name = name.replace('_', ' ')

        print(f'MMCD-{molecule_id}, {spectrum_id}, {self.shifts_as_list(shifts)}', file=self.id_shifts_csv_file)
        print(f'MMCD-{molecule_id}, "{name}"', file=self.id_name_csv_file)


if __name__ == "__main__":
    directory = '/home/mh491/Metameta_Files/hmdb_nmr_spectra'
    metabolite_reader = HMDB_Metabolites_to_CSV(directory)
    metabolite_reader.run()

    # directory = '/home/mh491/Metameta_Files/hmdb_nmr_spectra'
    # reader = HMDB_Reader(directory)
    # reader = HMDB_to_MYSQL(directory)
    # reader.run()

    # directory = '/home/mh491/Metameta_Files/bmrb_nmr_spectra'
    # reader = BMRB_to_CSV(directory)
    # reader.run()

    # directory = '/home/mh491/Metameta_Files/mmcd_nmr_spectra'
    # reader = MMCD_Reader(directory)
    # reader = MMCD_to_CSV(directory)
    # reader.run()
