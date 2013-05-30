import csv, logging

from collections import Counter

from pbhla.utils import BlasrM1
from pbhla.io.SamIO import SamReader

MIN_COUNT = 100000
MIN_FRAC = 0.05

class ReferenceSelector( object ): 

    def __init__(self, blasr_files, min_count=None, min_frac=None):
        self.blasr_files = blasr_files
        self.min_count = min_count if min_count else MIN_COUNT
        self.min_frac = min_frac if min_frac else MIN_FRAC
        self.selected_references = []
        self.initialize_logger()
        self.run()

    def run(self):
        for blasr_file in self.blasr_files:
            if blasr_file.endswith('.m1'):
                self.parse_blasr_file( blasr_file )
            else:
                msg = 'Unrecognized Blasr file-type! "{0}"'.format( blasr_file )
                self.log.info( msg )
                raise ValueError( msg )

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusReferenceSelector with the following:")
        self.log.info("\tMinimum Count: {0}".format(self.min_count))
        self.log.info("\tMinimum Fraction: {0}".format(self.min_frac))
        for filename in self.blasr_files:
            self.log.info("\tInput File: {0}".format(filename))

    def parse_blasr_file(self, blasr_file):
        reference_hits = self.parse_reference_hits( blasr_file )
        best_hits = self.identify_best_hits( reference_hits )
        self.selected_references +=  best_hits

    def parse_reference_hits(self, blasr_file):
        msg = 'Reading Blasr results from "{0}"'.format(blasr_file)
        self.log.info( msg )
        reference_hits = []
        with open(blasr_file, 'r') as handle:
            for hit in map(BlasrM1._make, csv.reader(handle, delimiter=' ')):
                reference_hits.append( hit.tname )
        self.log.info('Finished reading Blasr results')
        return reference_hits

    def identify_best_hits(self, reference_hits):
        self.log.info('Identifying best reference sequences')
        counts = Counter( reference_hits )
        top_two = counts.most_common(2)
        if len(top_two) == 0:
            return []
        if len(top_two) == 1:
            return [top_two[0][0]]
        first_ref, first_count = top_two[0]
        second_ref, second_count = top_two[1]
        if second_count < self.min_count:
            return [first_ref]
        elif second_count/float(first_count) < self.min_frac:
            return [first_ref]
        else:
            return [first_ref, second_ref]
        self.log.info('Finished identifying the best references')

    def __setitem__(self, key, value):
        self.selected_references[key] = value

    def __getitem__(self, key):
        return self.selected_references[key]

    def __iter__(self):
        return iter(self.selected_references)

    def __len__(self):
        return len(self.selected_references)
