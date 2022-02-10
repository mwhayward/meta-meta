import pathlib
import xml.etree.ElementTree as et
import nmrstarlib
import sys


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
        print(f'HMDB-{self.id_to_index(id)}, "{name}"')

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
                    name = name.text
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
        print(f'{id}, "{name}"', file=self.csv_file)


class HMDB_Reader:

    def __init__(self, directory):
        self._directory = directory

    def shifts_as_list(self, peak_list):
        return ', '.join(["%7.3f" % peak[0] for peak in peak_list])

    def intensities_as_list(self, peak_list):
        return ', '.join(["%7.3f" % peak[-1] for peak in peak_list])

    def output(self, file, molecule_id, spectrum_id, frequency, ph, reference, temperature, peak_list):
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
                    frequency = spectrum.find('frequency').text
                    ph = spectrum.find('sample-ph').text
                    reference = spectrum.find('chemical-shift-reference').text
                    temperature = spectrum.find('sample-temperature').text
                    peak_list = []

                    peaks = spectrum.find("nmr-one-d-peaks")

                    for peak in peaks:
                        shift = float(peak.find("chemical-shift").text)
                        intensity = float(peak.find("chemical-shift").text)
                        peak_list.append((shift, intensity))

                    peak_list = sorted(peak_list)

                    self.output(file, molecule_id, spectrum_id, frequency, ph, reference, temperature, peak_list)


class HMDB_to_CSV(HMDB_Reader):

    def __init__(self, directory):
        super(HMDB_to_CSV, self).__init__(directory)
        self.metabolites_out_file = 'hmdb_xml_metbolites.csv'
        self. multiplets_out_file = 'hmdb_xml_mmultiplets.csv'
        self.peaks_out_file = 'hmdb_xml_peaks.csv'

    def run(self):
        with open(pathlib.Path(directory, self.metabolites_out_file), 'w') as self.metabolites_csv_file, \
             open(pathlib.Path(directory, self.multiplets_out_file), 'w') as self.multiplets_csv_file, \
             open(pathlib.Path(directory, self.peaks_out_file), 'w') as self.peaks_csv_file:
            super(HMDB_to_CSV, self).run()

    def output(self, file, molecule_id, spectrum_id, frequency, ph, reference, temperature, peak_list):
        shifts_text = self.shifts_as_list(peak_list)
        intensities_text = self.intensities_as_list(peak_list)
        print(f'{molecule_id}, {spectrum_id}, {frequency}, {ph}, {reference}, {temperature}', file=self.metabolites_csv_file)
        for i, peak in enumerate(peak_list):
            print(f'{molecule_id}.{i+1}, {molecule_id}, {peak[0]}, {peak[-1]}', file=self.peaks_csv_file)

        # TODO: complete output to write to each file for spectrum, multiplets and peaks.
        #  Might not need the shifts/intensities_text methods


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
                                if frame[key] == 'theoretical chemical shifts' and result != True:
                                    result = True
        return result

    def get_1h_shifts(self, tree):
        result = []
        SAVE_FRAME = 'save_theoretical_chem_shifts'
        LOOP_FRAME = 'Theoretical_chem_shift.Atom_type'
        for top_frame in [top for top in tree if top.startswith(SAVE_FRAME)]:
            for loop in [elem for elem in tree[top_frame] if elem.startswith('loop')]:
                for elem in tree[top_frame][loop][1:]:
                    for i, line in enumerate(elem):
                        for key in line.keys():
                            if key == LOOP_FRAME and line[LOOP_FRAME] == 'H':
                                try:
                                    result.append(float(line['Theoretical_chem_shift.Val']))
                                except:
                                    print(f"   cound\'t convert {line['Theoretical_chem_shift.Val']} to float")
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
            file = pathlib.Path(directory, f'bmst{str(id).zfill(6)}.str')

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
    # directory = '/home/mh491/Metameta_Files/hmdb_nmr_spectra'
    # metabolite_reader = HMDB_Metabolites_to_CSV(directory)
    # metabolite_reader.run()

    directory = '/home/mh491/Metameta_Files/hmdb_nmr_spectra'
    # reader = HMDB_Reader(directory)
    reader = HMDB_to_CSV(directory)
    reader.run()

    # directory = '/home/mh491/Metameta_Files/bmrb_nmr_spectra'
    # reader = BMRB_to_CSV(directory)
    # reader.run()

    # directory = '/home/mh491/Metameta_Files/mmcd_nmr_spectra'
    # reader = MMCD_Reader(directory)
    # reader = MMCD_to_CSV(directory)
    # reader.run()
