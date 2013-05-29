import re
import csv
import logging

from pbcore.io.FastaIO import FastaReader, FastaWriter

class LocusReference( object ):
    """
    A class for tabulating subread statistics for a single HLA Locus
    """

    def __init__(self, input_file, output_file, header=False, id_list=None):
        self.input_file = input_file
        self.output_file = output_file
        self.header = header
        self.id_list = set(id_list) if id_list else None
        self.initialize_logger()
        self.create_locus_reference()

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusReference with the following:")
        self.log.info("\tInput File: {0}".format(self.input_file))
        self.log.info("\tOutput File: {0}".format(self.output_file))
        self.log.info("\tHeader: {0}".format(self.header))
        self.log.info("\tId List: {0}".format(self.id_list))

    def create_locus_reference(self):
        cumulative_records = []
        with open(self.input_file, 'r') as handle:
            for line in handle:
                filename, locus = line.strip().split()
                records = list( FastaReader( filename ) )
                if self.id_list:
                    records = [r for r in records if r.name.split()[0] in self.id_list]
                if self.header:
                    records = [records[0]]
                cumulative_records += records
        self.output_records( cumulative_records )

    def output_records(self, records):
        with FastaWriter( self.output_file ) as handle:
            for record in records:
                handle.writeRecord( record )
