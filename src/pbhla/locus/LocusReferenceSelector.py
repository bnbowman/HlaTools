import csv, logging

from collections import Counter

from pbhla.utils import BlasrM1

MIN_COUNT = 100
MIN_FRAC = 0.05

class LocusReferenceSelector( object ): 

    def __init__(self, blasr_file, min_count=None, min_frac=None):
        self.blasr_file = blasr_file
        self.min_count = min_count if min_count else MIN_COUNT
        self.min_frac = min_frac if min_frac else MIN_FRAC
        self.reference_list = []
        self.initialize_logger()
        if self.blasr_file.endswith('.m1'):
            self.parse_blasr_file()
        else:
            msg = 'Unrecognized Blasr file-type! "{0}"'.format(self.blasr_file)
            self.log.info( msg )
            raise ValueError( msg )
        return self.identify_best_references()

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusDict with the following:")
        self.log.info("\tBlasr File: {0}".format(self.blasr_file))
        self.log.info("\tMinimum Count: {0}".format(self.min_count))
        self.log.info("\tMinimum Fraction: {0}".format(self.min_frac))

    def parse_blasr_file(self):
        msg = 'Reading Blasr results from "{0}"'.format(self.blasr_file)
        self.log.info( msg )
        with open(self.blasr_file, 'r') as handle:
            for hit in map(BlasrM1._make, csv.reader(handle, delimiter=' ')):
                self.reference_hits.append( hit.tname )
        self.log.info('Finished reading Blasr results')

    def identify_best_references(self):
        self.log.info('Identifying best reference sequences')
        counts = Counter(self.reference_list)
        top_two = counts.most_common(2)
        first_ref, first_count = top_two[0]
        second_rec, second_count = top_two[1]
        if second_count < self.min_count:
            return [first_ref]
        elif second_count/float(first_count) < self.min_frac:
            return [first_ref]
        else:
            return [first_ref, second_ref]
