# /usr/bin/env python

import os, logging

from pbhla.io.BasH5Reader import BasH5Reader

MIN_LENGTH = 500
MIN_SCORE = 0.75

class BasH5Extractor( object ):
    """
    A tool for extracting sub-reads from BasH5 files
    """

    def __init__( self, input_file,
                        output_file,
                        min_length=None,
                        min_score=None,
                        dilution=None ):
        self.input_file = input_file
        self.output_file = output_file
        self.min_length = min_length if min_length else MIN_LENGTH
        self.min_score = min_score if min_score else MIN_SCORE
        self.dilution = dilution
        self.initialize_logger()
        self.run()

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info('Extracting reads according to the following criteria:')
        self.log.info('\t\tMinimum Subread Length: %s' % self.min_length)
        self.log.info('\t\tMinimum Read Score: %s' % self.min_score)
        self.log.info('\t\tSubread Dilution Factor: %s' % self.dilution)

    def run(self):
        file_list = self.parse_input_file()
        self.extract_all_subreads( file_list )

    def check_filepath(self, filepath):
        if not os.path.isfile( filepath ):
            msg = 'Invalid input file for subread extraction "{0}"!'.format(filepath)
            self.log.info( msg )
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
