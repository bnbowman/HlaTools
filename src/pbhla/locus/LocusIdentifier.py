import csv, logging

from pbhla.utils import BlasrM1

class LocusIdentifier( object ): 

    def __init__(self, blasr_file, reference_dict):
        self.blasr_file = blasr_file
        self.reference_dict = reference_dict
        self.locus_dict = {}
        self.initialize_logger()
        if self.blasr_file.endswith('.m1'):
            self.parse_blasr_file()
        else:
            msg = 'Unrecognized Blasr file-type! "{0}"'.format(self.blasr_file)
            self.log.info( msg )
            raise ValueError( msg )

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusDict with the following:")
        self.log.info("\tBlasr File: {0}".format(self.blasr_file))

    def parse_blasr_file(self):
        msg = 'Reading Blasr results from "{0}"'.format(self.blasr_file)
        self.log.info( msg )
        with open(self.blasr_file, 'r') as handle:
            for hit in map(BlasrM1._make, csv.reader(handle, delimiter=' ')):
                locus = self.reference_dict[hit.tname]
                self.add_record( hit.qname, locus )
        self.log.info('Finished reading Blasr results')

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
