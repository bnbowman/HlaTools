import logging

from pbcore.io.FastaIO import FastaReader

class LocusDict( object ): 

    def __init__(self, input_file):
        self.input_file = input_file
        self.locus_dict = {}
        self.initialize_logger()
        if self.input_file.endswith('.fofn'):
            self.parse_reference_fofn()
        elif self.input_file.endswith('.txt'):
            self.parse_locus_key()
        else:
            msg = 'Unrecognized input file-type! "{0}"'.format(self.input_file)
            self.log.info( msg )
            raise ValueError( msg )

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusDict with the following:")
        self.log.info("\tInput File: {0}".format(self.input_file))

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
                    self.add_record( name, locus )
        self.log.info('Finished reading Locus References')

    def parse_locus_key(self):
        msg = 'Reading existing Locus Key from "{0}"'.format(self.input_file)
        self.log.info( msg )
        with open(self.input_file, 'r') as handle:
            for line in handle:
                seq_name, locus = line.strip().split()
                self.add_record( seq_name, locus )
        self.log.info('Finished reading Locus Key')

    def add_record(self, name, locus):
        if name in self.locus_dict:
            msg = 'Duplicate sequence ids found! "{0}"'.format( name )
            self.log.info( msg )
            raise KeyError( msg )
        self.locus_dict[name] = locus

    def write(self, output_file):
        msg = 'Writing new Locus Key to "{0}"'.format(output_file)
        self.log.info( msg )
        with open(output_file, 'w') as handle:
            for seq_name, locus in self.locus_dict.iteritems():
                print >> handle, '{0} {1}'.format(seq_name, locus)

    def __setitem__(self, key, value):
        self.locus_dict[key] = value

    def __getitem__(self, key):
        return self.locus_dict[key]

    def __delitem__(self, key):
        del self.locus_dict[key]

    def __iter__(self):
        return iter(self.locus_dict)
