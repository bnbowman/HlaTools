import csv, logging

from pbhla.utils import BlasrM1
from pbhla.io.SamIO import SamReader
from pbcore.io.FastaIO import FastaReader

class ReferenceDict( object ): 

    def __init__(self, input_file, reference=None):
        self._dict = {}
        self.input_file = input_file
        self.file_type = self.input_file.split('.')[-1]
        self.reference = reference
        self.initialize_logger()
        self.validate_input()
        self.run()

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusDict with the following:")
        self.log.info("\tInput File: {0}".format(self.input_file))
        self.log.info("\tReference: {0}".format(self.reference))

    def validate_input(self):
        if self.file_type in ['fofn', 'txt', 'm1', 'sam']:
            pass
        else:
            msg = 'Unrecognized Alignment file-type! "{0}"'.format(self.input_file)
            self.log.info( msg )
            raise ValueError( msg )

    def run(self):
        if self.file_type == 'fofn':
            self.parse_reference_fofn()
        if self.file_type == 'txt':
            self.parse_locus_key()
        if self.file_type == 'm1':
            self.parse_blasr_file()
        if self.file_type == 'sam':
            self.parse_sam_file()

    def parse_reference_fofn(self):
        msg = 'Reading Locus References from "{0}"'.format(self.input_file)
        self.log.info( msg )
        with open(self.input_file, 'r') as handle:
            for line in handle:
                fasta_file, locus = line.strip().split()
                msg = 'Reading "{0}" sequences from "{1}"'.format(locus, fasta_file)
                self.log.info( msg )
                for record in FastaReader( fasta_file ):
                    name = record.name.split()[0]
                    name = name.split('_')[0]
                    self[name] = locus
        self.log.info('Finished reading Locus References')

    def parse_locus_key(self):
        msg = 'Reading existing Locus Key from "{0}"'.format(self.input_file)
        self.log.info( msg )
        with open(self.input_file, 'r') as handle:
            for line in handle:
                seq_name, locus = line.strip().split()
                self[seq_name] = locus
        self.log.info('Finished reading Locus Key')

    def parse_blasr_file(self):
        msg = 'Parsing Blasr results from "{0}"'.format(self.input_file)
        self.log.info( msg )
        with open(self.input_file, 'r') as handle:
            for hit in map(BlasrM1._make, csv.reader(handle, delimiter=' ')):
                key = hit.qname.split('/')[0]
                if self.reference:
                    self[key] = self.reference[hit.tname]
                else:
                    self[key] = hit.tname
        self.log.info('Finished reading Blasr results')

    def parse_sam_file(self):
        msg = 'Parsing SAM alignments from "{0}"'.format(self.input_file)
        self.log.info( msg )
        for record in SamReader(self.input_file):
            key = record.qname.split('/')[0]
            if self.reference:
                self[key] = self.reference[record.rname]
            else:
                self[key] = record.rname
        self.log.info('Finished reading SAM file results')

    def write(self, output_file):
        msg = 'Writing new Locus Key to "{0}"'.format(output_file)
        self.log.info( msg )
        with open(output_file, 'w') as handle:
            for seq_name, locus in self._dict.iteritems():
                print >> handle, '{0} {1}'.format(seq_name, locus)

    def __setitem__(self, key, value):
        if key in self._dict:
            msg = 'Duplicate sequence ids found! "{0}"'.format( key )
            self.log.info( msg )
            raise KeyError( msg )
        self._dict[key] = value

    def __getitem__(self, key):
        return self._dict[key]

    def __delitem__(self, key):
        del self._dict[key]

    def __iter__(self):
        return iter(self._dict)
