import os, logging

from pbcore.io.FastaIO import FastaReader, FastaWriter

class SequenceSeparator( object ): 

    def __init__(self, fasta_file, reference_dict=None, selected=None):
        self.fasta_file = fasta_file
        self.reference_dict = reference_dict
        self.selected = selected
        self.sequences = {}
        self.initialize_logger()
        self.validate_input_file()
        self.run()

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing LocusSeparator with the following:")
        self.log.info("\tFasta File: {0}".format(self.fasta_file))
        self.log.info("\tSelected Values: {0}".format(self.selected))

    def validate_input_file( self ):
        if self.fasta_file.endswith('.fa') or self.fasta_file.endswith('.fasta'):
            pass
        else:
            msg = 'Unrecognized Fasta file-type! "{0}"'.format(self.fasta_file)
            self.log.info( msg )
            raise ValueError( msg )

    def run(self):
        if self.reference_dict and self.selected:
            self.separate_by_selected_reference()
        elif self.selected:
            self.separate_by_selected_id()
        elif self.reference_dict:
            self.separate_by_reference()

    def separate_by_selected_reference(self):
        pass

    def separate_by_selected_id(self):
        self.log.info('Separating selected references')
        for record in FastaReader( self.fasta_file ):
            name = record.name.split()[0]
            if name in self.selected:
                self.add_record( record, 'selected' )
            else:
                self.add_record( record, 'not_selected' )
        self.log.info('Finished separating sequences')

    def separate_by_reference(self):
        self.log.info('Separating sequences by reference')
        for record in FastaReader( self.fasta_file ):
            name = record.name.split()[0]
            ref = self.get_reference( name )
            self.add_record( record, ref )
        self.log.info('Finished separating sequences')

    def get_reference(self, name):
        try:
            reference = self.reference_dict[name]
        except:
            reference = 'Unmapped'
        return reference

    def add_record(self, record, locus):
        try:
            self.sequences[locus].append( record )
        except:
            self.sequences[locus] = [ record ]

    def write(self, value, output_file):
        self.log.info('Writing "{0}" sequences to "{1}"'.format(value, output_file))
        with FastaWriter( output_file ) as handle:
            for record in self.sequences[value]:
                handle.writeRecord( record )
        self.log.info('Finished writing sequences from "{0}"'.format(value))
        return output_file

    def write_all(self, prefix='Locus'):
        self.log.info('Writing fasta sequences from all loci to file')
        output_files = []
        for value in self.sequences:
            output_file = '{0}_{1}.fasta'.format(prefix, value)
            output_file = self.write(value, output_file)
            output_files.append( output_file )
        self.log.info('Finished writing sequences from all loci')
        return output_files
