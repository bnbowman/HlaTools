# /usr/bin/env python

import os, logging

from pbcore.io.BasH5Reader import BasH5Reader

MIN_LENGTH = 1000
MIN_SCORE = 0.80
DILUTION = 1.0

class SubreadExtractor( object ):
    """
    A tool for extracting sub-reads from BasH5 files
    """

    def __init__( self, input_file,
                        output_file,
                        min_length=None,
                        min_score=None,
                        dilution=None ):
        # Set input and output files
        self.input_file = input_file
        self.output = os.path.abspath( output )
        self.output_fasta = self.output + '.fasta'
        self.output_text = self.output + '.txt'
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
        file_list = self.parse_input_file()
        self.extract_all_subreads( file_list )

    def check_filepath(self, filepath):
        if not os.path.isfile( filepath ):
            msg = 'Invalid input file for subread extraction "{0}"!'.format(filepath)
            log.info( msg )
            raise ValueError( msg )

    def parse_input_file( self ):
        """Convert the input file to a list of absolute paths to Bash5s"""
        self.check_filepath( self.input_file )
        if self.input_file.endswith('.bas.h5'):
            return [ self.input_file ]
        elif self.input_file.endswith('.fofn'):
            file_list = []
            with open( self.input_file, 'r' ) as handle:
                for line in handle:
                    filepath = line.strip()
                    self.check_filepath( filepath )
                    file_list.append( filepath )
            return file_list
        else:
            msg = "Input file must be Bas.H5 of FOFN!"
            self.log.info( msg )
            raise TypeError( msg )

    def extract_all_subreads( self, file_list ):
        with open(self.output_file, 'w') as handle:
            for filename in file_list:
                self.write_subreads( filename, handle )

    def write_subreads(self, filename, handle):
        for zmw in BasH5Reader( filename ):
            if zmw.readScore < self.min_score:
                continue
            adapter_starts = [z.readStart for z in zmw.adapters]
            adapter_ends = [z.readEnd for z in zmw.adapters]
            for subread in zmw.subreads:
                if len(subread.basecalls()) < self.min_length:
                    continue
                if subread.readStart in adapter_ends and subread.readEnd in adapter_starts:
                    seq_name = "fp_{0}_{1}_{2}".format(subread.holeNumber, subread.readStart, subread.readEnd)
                else:
                    seq_name = "{0}_{1}_{2}".format(subread.holeNumber, subread.readStart, subread.readEnd)
                print >> handle, ">{0}\n{1}".format(seq_name, subread.basecalls())

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
