from os import listdir
from os.path import isfile, join
import pandas as pd

metabolites = {'accessions': [],
               'numbers':    []}

path = '/home/mh491/Database/HMDB_files/nmrML_experimental_Feb15_2022'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
for file in onlyfiles:
    if file.split('_')[2] == '1H.nmrML':
        acc = file.split('_')[0]
        number = file.split('_')[1]
        if acc not in metabolites['accessions']:
            metabolites['accessions'].append(acc)
            metabolites['numbers'].append(number)

path = '/home/mh491/Database/HMDB_files/hmdb_nmr_peak_lists'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
for file in onlyfiles:
    if file.split('_')[1] == 'nmroned':
        acc = file.split('_')[0]
        number = file.split('_')[2]
        if acc not in metabolites['accessions']:
            metabolites['accessions'].append(acc)
            metabolites['numbers'].append(number)

path = '/home/mh491/Database/HMDB_files/xml_files'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
for file in onlyfiles:
    if file.endswith('xml') and file.split('_')[2] == 'one':
        acc = file.split('_')[0]
        number = file.split('_')[5].rstrip('.xml')
        if acc not in metabolites['accessions']:
            metabolites['accessions'].append(acc)
            metabolites['numbers'].append(number)

metabolites = pd.DataFrame(metabolites['numbers'], index=metabolites['accessions'], columns=['number'])

id_names = pd.read_csv('/home/mh491/Database/HMDB_files/xml_files/id_name.csv', index_col=0, header=None, names=['Name'])

merged = metabolites.join(id_names)
final = merged.drop_duplicates(subset=['Name'])
final.to_csv('/home/mh491/Database/HMDB_files/metabolite_names.csv')



