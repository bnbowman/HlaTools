#! /usr/bin/env python

import os, logging

from pbcore.io.BasH5Reader import BasH5Reader
from pbcore.io.FastaIO import FastaRecord, FastaWriter

MIN_LENGTH = 2000
MIN_SCORE = 0.80
DILUTION = 1.0

log = logging.getLogger()

class SubreadExtractor( object ):
    """
    A tool for extracting sub-reads from BasH5 files
    """

    def __init__( self, input_file,
                        output,
                        min_length=None,
                        min_score=None,
                        dilution=None ):
        # Set input and output files
        self.input_file = input_file
        self.output = os.path.abspath( output )
        self.fasta_file = self.output + '.fasta'
        self.text_file = self.output + '_fullpass.txt'
        # Set other variables
        self.min_length = min_length or MIN_LENGTH
        self.min_score = min_score or MIN_SCORE
        self.dilution = dilution or DILUTION
        # Log the settings
        log.info('Extracting reads according to the following criteria:')
        log.info('\t\tMinimum Length: %s' % self.min_length)
        log.info('\t\tMinimum Score: %s' % self.min_score)
        log.info('\t\tDilution Factor: %s' % self.dilution)

    def __call__(self):
        # Validate and parse the input file
        check_filepath( self.input_file )
        file_list = parse_input_file( self.input_file )
        # Extract subreads from each of the identified file
        with FastaWriter( self.fasta_file ) as fasta_handle:
            with open( self.text_file, 'w' ) as text_handle:
                for filename in file_list:
                    self.extract_subreads( filename, fasta_handle, text_handle )
        # Return the paths to the output files
        return self.fasta_file, self.text_file

    def extract_subreads(self, filename, fasta_handle, text_handle):
        """
        Extract sub-reads from a file and 
        """
        for zmw in BasH5Reader( filename ):
            if zmw.readScore < self.min_score:
                continue
            adapter_starts = [z.readStart for z in zmw.adapters]
            adapter_ends = [z.readEnd for z in zmw.adapters]
            for subread in zmw.subreads:
                if len(subread.basecalls()) < self.min_length:
                    continue
                record = FastaRecord( subread.readName, subread.basecalls() )
                fasta_handle.writeRecord( record )
                # If the read goes from adaptor to adaptor, record it as FP
                if ((subread.readStart in adapter_ends) and
                    (subread.readEnd in adapter_starts)):
                    text_handle.write( record.name + '\n' )

def parse_input_file( input_file ):
    """
    Convert the input file to a list of absolute paths to Bash5s
    """
    if input_file.endswith('.bas.h5') or input_file.endswith('bax.h5'):
        return [ input_file ]
    elif input_file.endswith('.fofn'):
        file_list = []
        with open( input_file, 'r' ) as handle:
            for line in handle:
                filepath = line.strip()
                check_filepath( filepath )
                file_list.append( filepath )
        return file_list
    else:
        msg = "Input file must be Bas.H5 of FOFN!"
        log.info( msg )
        raise TypeError( msg )

def check_filepath( filepath ):
    """
    Check that a specified file exists
    """
    if not os.path.isfile( filepath ):
        msg = 'Invalid input file for subread extraction "{0}"!'.format(filepath)
        log.info( msg )
        raise ValueError( msg )

if __name__ == '__main__':
    import argparse
    desc = 'Extract subreads and identify fullpass reads'
    parser = argparse.ArgumentParser( description=desc )

    add = parser.add_argument
    add('input_file', metavar='INPUT',
        help='Input Bas.H5 or FOFN file to extract reads from')
    add('output', metavar='PREFIX',
        help='Output prefix for both reads (Fasta) and Pass-Status (Txt)')
    add('--min_length', type=int, default=MIN_LENGTH,
        help='Minimum length for subreads ({0})'.format(MIN_LENGTH))
    add('--min_score', type=float, default=MIN_SCORE,
        help='Minimum ReadScore for subreads ({0})'.format(MIN_SCORE))
    add('--dilution', type=float, default=DILUTION,
        help='Fraction of reads to use ({0})'.format(DILUTION))
    args = parser.parse_args()

    extractor = SubreadExtractor( args.input_file,
                                  args.output,
                                  args.min_length,
                                  args.min_score,
                                  args.dilution )
    extractor()
