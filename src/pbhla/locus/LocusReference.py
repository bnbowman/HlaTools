#!/home/UNIXHOME/jquinn/HGAP_env/bin/python
import re
import csv
import logging

from pbhla.io.FastaIO import FastaReader, FastaWriter

class LocusReference( object ):
    """
    A class for tabulating subread statistics for a single HLA Locus
    """

    def __init__(self, input_file, output_file, header=False):
        self.input_file = input_file
        self.output_file = output_file
        self.header = header
        self.initialize_logger()
        self.create_locus_reference()

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusReference with the following:")
        self.log.info("\tInput File: {0}".format(self.input_file))
        self.log.info("\tOutput File: {0}".format(self.output_file))

    def create_locus_reference(self):
        sequences = {}
        with open(self.input_file, 'r') as handle:
            for line in handle:
                filename, locus = line.strip().split()
                records = list( FastaReader( filename ) )
        self.output_records( records )

    def output_records(self, records):
        with FastaWriter( self.output_file ) as handle:
            for record in records:
                handle.writeRecord( record )
