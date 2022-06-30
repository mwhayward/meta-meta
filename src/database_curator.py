from pathlib import Path
import statistics
import xml.etree.ElementTree as et
import sys
import os
import pandas as pd
import sqlite3
from sqlite3 import Error
import re


class Curator:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.conn = self.create_connection()

    def create_connection(self):
        # Creates a single connection to the database
        # This connection is used at the end of the script to export the tables to a db file
        conn = None
        try:
            conn = sqlite3.connect(str(self.directory.joinpath('ccpn_metabolites.db')))
        except Error as e:
            print(e)
        return conn

    def run(self):
        tables = {}
        for table in ['metabolites', 'synonyms', 'isin', 'ontology', 'concentrations',
                      'samples', 'spectra', 'multiplets', 'peaks']:
            query = f'select * from {table}'
            tables[table] = pd.read_sql(query, self.conn)
        for table in tables:
            print(table)
            print(tables[table])


if __name__ == "__main__":
    directory = '/home/mh491/Database'
    curator = Curator(directory)
    curator.run()
