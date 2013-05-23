import os, logging

from pbhla.io.FastaIO import FastaReader, FastaWriter

class SubreadSeparator( object ): 

    def __init__(self, fasta_file, reference_dict):
        self.fasta_file = fasta_file
        self.reference_dict = reference_dict
        self.locus_sequences = {}
        self.initialize_logger()
        if self.fasta_file.endswith('.fa') or self.fasta_file.endswith('.fasta'):
            self.separate_fasta_file()
        else:
            msg = 'Unrecognized Fasta file-type! "{0}"'.format(self.fasta_file)
            self.log.info( msg )
            raise ValueError( msg )

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusSeparator with the following:")
        self.log.info("\tFasta File: {0}".format(self.fasta_file))

    def separate_fasta_file(self):
        msg = 'Separating sequences from "{0}"'.format(self.fasta_file)
        self.log.info( msg )
        for record in FastaReader( self.fasta_file ):
            locus = self.get_sequence_locus( record.name )
            self.add_record( record, locus )
        msg = 'Finished separating sequences'
        self.log.info( msg )

    def get_sequence_locus(self, name):
        try:
            locus = self.reference_dict[name]
        except:
            locus = 'Unmapped'
        return locus

    def add_record(self, record, locus):
        try:
            self.locus_sequences[locus].append( record )
        except:
            self.locus_sequences[locus] = [ record ]

    def write(self, locus, prefix='Locus'):
        output_file = '{0}_{1}.fasta'.format(prefix, locus)
        if os.path.isfile( output_file ):
            self.log.info('Found existing output file "{0}"'.format(output_file))
            self.log.info('Skipping alignment step for "{0}"'.format(locus))
            return output_file
        msg = 'Writing "{0}" sequences to "{1}"'.format(locus, output_file)
        self.log.info( msg )
        with FastaWriter( output_file ) as handle:
            for record in self.locus_sequences[locus]:
                handle.writeRecord( record )
        self.log.info('Finished writing sequences from "{0}"'.format(locus))
        return output_file

    def write_all(self, prefix='Locus'):
        self.log.info('Writing fasta sequences from all loci to file')
        output_files = []
        for locus in self.locus_sequences:
            output_file = self.write(locus, prefix)
            output_files.append( output_file )
        self.log.info('Finished writing sequences from all loci')
        return output_files
