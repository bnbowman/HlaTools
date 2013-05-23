import csv, logging

from pbhla.utils import BlasrM1
from pbhla.io.SamIO import SamReader

class SubreadLocusDict( object ): 

    def __init__(self, input_file, reference_dict):
        self.input_file = input_file
        self.reference_dict = reference_dict
        self.locus_dict = {}
        self.initialize_logger()
        self.run()

    def run(self):
        if self.input_file.endswith('.m1'):
            self.parse_blasr_file()
        elif self.input_file.endswith('.sam'):
            self.parse_sam_file()
        else:
            msg = 'Unrecognized Alignment file-type! "{0}"'.format(self.input_file)
            self.log.info( msg )
            raise ValueError( msg )

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusDict with the following:")
        self.log.info("\tInput File: {0}".format(self.input_file))

    def parse_blasr_file(self):
        msg = 'Parsing Blasr results from "{0}"'.format(self.input_file)
        self.log.info( msg )
        with open(self.input_file, 'r') as handle:
            for hit in map(BlasrM1._make, csv.reader(handle, delimiter=' ')):
                locus = self.reference_dict[hit.tname]
                query = hit.qname.split('/')[0]
                self.add_record( hit.qname, locus )
        self.log.info('Finished reading Blasr results')

    def parse_sam_file(self):
        msg = 'Parsing SAM alignments from "{0}"'.format(self.input_file)
        self.log.info( msg )
        for record in SamReader(self.input_file):
            locus = self.reference_dict[record.rname]
            query = record.qname.split('/')[0]
            self.add_record( query, locus )
        self.log.info('Finished reading SAM file results')

    def add_record(self, name, locus):
        if name in self.locus_dict:
            msg = 'Duplicate sequence ids found! "{0}"'.format( name )
            self.log.info( msg )
            raise KeyError( msg )
        self.locus_dict[name] = locus

    def __setitem__(self, key, value):
        self.locus_dict[key] = value

    def __getitem__(self, key):
        return self.locus_dict[key]

    def __delitem__(self, key):
        del self.locus_dict[key]

    def __iter__(self):
        return iter(self.locus_dict)
