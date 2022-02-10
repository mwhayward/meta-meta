import requests
import urllib.request
import zipfile
import pathlib
from html.parser import HTMLParser




class HMDB_Downloader:

    def get_file_name_from_url(self, url):
        url_path = pathlib.PurePosixPath(urllib.parse.urlparse(url).path)
        return url_path

    def run(self, directory, biofluid='hmdb'):

        URL_SPECTRA = 'http://specdb.wishartlab.com/downloads/exports/spectra_xml/hmdb_nmr_spectra.zip'
        URL_METABOLITES = 'http://www.hmdb.ca/system/downloads/current/' + biofluid + '_metabolites.zip'

        spectra_name =self.get_file_name_from_url(URL_SPECTRA).stem
        target_directory = pathlib.Path(directory, spectra_name)

        self.download_and_extract_xml_zip(URL_SPECTRA, target_directory)
        self.download_and_extract_xml_zip(URL_METABOLITES, target_directory)



    def download_and_extract_xml_zip(self, URL, directory):

        file = self.get_file_name_from_url(URL).name
        zip_target = pathlib.Path(directory, file)

        target = pathlib.Path(directory)
        if not target.exists():
            target.mkdir()
        print(f'Beginning download of {file}...')
        print(f"    target directory is {target}")

        urllib.request.urlretrieve(URL, zip_target)
        print(f'{file} successfully downloaded')

        zip_file = pathlib.Path(directory, file)
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(target)

        print(f"download and extract complete")
        print(f"    removing zip file  {zip_target}")

        zip_target.unlink()


class BMRB_Downloader:
    class BMRB_Directory_Parser(HTMLParser):

        def __init__(self):
            super().__init__()
            self.files = []

        def handle_starttag(self, tag, attrs):
            if tag == 'a':
                href = [elem for elem in attrs if elem[0] == 'href'][0][1]
                if href.startswith('bmse'):
                    file_name = f'{attrs[0][1][:-1]}.str'
                    href = f'{href}{file_name}'
                    self.files.append(href)

    def run(self,directory):
        URL = 'https://bmrb.io/ftp/pub/bmrb/metabolomics/entry_directories/'
        bmrb_directory = urllib.request.urlopen(URL).read().decode('utf-8')

        parser = BMRB_Downloader.BMRB_Directory_Parser()
        parser.feed(bmrb_directory)

        download_directory = pathlib.Path(directory, 'bmrb_nmr_spectra')
        download_directory.mkdir(exist_ok=True)

        print('xml download start:')

        for i, file in enumerate(parser.files):
            file_url = URL + file
            bruker_file_url = URL + pathlib.Path(file).parts[0] + '/nmr'
            star_file_end = pathlib.Path(file).parts[-1]
            xml_file_end = pathlib.Path(file).parts[0] + '.xml'
            star_target = pathlib.Path(download_directory, star_file_end)
            xml_target = pathlib.Path(download_directory, xml_file_end)
            # urllib.request.urlretrieve(file_url, target)
            response = requests.get(bruker_file_url +'/set01')
            if response.status_code == 200:
                response = requests.get(bruker_file_url +'/set01/1/pdata/1/peaklist.xml')
                if response.status_code == 200:
                    urllib.request.urlretrieve(bruker_file_url + '/set01/1/pdata/1/peaklist.xml', xml_target)
                else:
                    response = requests.get(bruker_file_url + '/set01/1H/pdata/1/peaklist.xml')
                    if response.status_code == 200:
                        urllib.request.urlretrieve(bruker_file_url + '/set01/1H/pdata/1/peaklist.xml', xml_target)
                    else:
                        print(f'no bruker file found for {xml_file_end}')

            else:
                print(f'no {xml_file_end}')
                '''try:
                    urllib.request.urlopen(bruker_file_url +'set01/1/pdata/1')
                    print(f'{file} found')
                except:
                    print(f'{bruker_file_url} not found')'''
            # print(f'download {file} {i+1} of {len(parser.files)}')


class MMCD_Downloader:
    # class BMRB_Directory_Parser(HTMLParser):
    #
    #     def __init__(self):
    #         super().__init__()
    #         self.files = []
    #
    #     def handle_starttag(self, tag, attrs):
    #         if tag == 'a':
    #             href  = [elem for elem in attrs if elem[0] == 'href'][0][1]
    #             if href.endswith('.str'):
    #                 self.files.append(href)

    def run(self,directory):
        URL = 'http://mmcd.nmrfam.wisc.edu/peaklist/'
        FILE_NAME_TEMPLATE = 'expnmr_%05i_3.txt'


        download_directory  = pathlib.Path(directory,'mmcd_nmr_spectra')
        download_directory.mkdir(exist_ok=True)

        for i in range(1,1000):
            file_name = FILE_NAME_TEMPLATE % i
            file_url = URL + file_name

            target = pathlib.Path(download_directory, file_name)
            print (f'download {file_name}')

            try:
                urllib.request.urlretrieve(file_url, target)
            except:
                print("  failed to download...")

if __name__ == '__main__':
    # downloader = HMDB_Downloader()
    # downloader.run('/home/mh491/Metameta_Files')

    downloader = BMRB_Downloader()
    downloader.run('/home/mh491/Metameta_Files')

    # downloader = MMCD_Downloader()
    # downloader.run('/home/mh491/Metameta_Files')