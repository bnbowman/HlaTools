#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
import logging, logging.config

from pbhla import __LOG__
from pbhla.resequencing.Resequencer import Resequencer
from pbhla.bash5 import get_bash5_reader, filter_zmw_list,  write_zmw_whitelist
from pbhla.barcodes import get_barcode_reader, get_barcodes, get_barcode_zmws
from pbhla.fofn import create_baxh5_fofn
from pbhla.sequences.input import get_input_file
from pbhla.io.AmpAnalysisIO import AmpliconAnalysisReader
from pbhla.io.utils import write_records, get_unique_records
from pbhla.utils import create_directory

logging.config.fileConfig( __LOG__ )
log = logging.getLogger( __name__ )

class AmpliconAnalysisResequencer( object ):
    """
    An object for iterating over all of the barcodes in
    """

    def __init__(self, output, setup=None, nproc=1, debug=False):
        """Initialize cross-cluster and object-specific settings"""
        if debug:
            log.setLevel( logging.DEBUG )
            log.debug("TESTING")
        log.info("Initializing Resequencer sub-module")
        self._resequencer = Resequencer(setup, nproc)

        # Initialize output folder
        self._output = output
        create_directory( self.output )

    @property
    def resequencer(self):
        return self._resequencer

    @property
    def setup(self):
        return self.resequencer.setup

    @property
    def nproc(self):
        return self.resequencer.nproc

    @property
    def output(self):
        return self._output

    def get_output_folder(self, barcode):
        output_dir = os.path.join( self.output, barcode )
        create_directory( output_dir )
        return output_dir

    def __call__(self, amp_analysis, data_file, barcode_file, barcode_string=None, min_snr=None, min_length=None):
        log.info("Beginning Amplicon Analysis resequencing workflow for {0}".format(amp_analysis))

        # Pick or create a single file from AA and read it
        amp_analysis_file = get_input_file( amp_analysis )
        amp_analysis_records = list(AmpliconAnalysisReader(amp_analysis_file))

        # Convert the raw data file into a BaxH5 fofn for use downstream
        # and create appropriate reader for local access
        bash5 = get_bash5_reader( data_file )
        baxh5_file = os.path.join( self.output, 'baxh5.fofn')
        create_baxh5_fofn( data_file, baxh5_file )

        # Create a Reader for the Barcode data and find the overlap with any
        # barcodes specified by the user
        bc_reader = get_barcode_reader( barcode_file )
        bc_list = get_barcodes( bc_reader, barcode_string )

        for i, bc in enumerate( bc_list ):
            log.info('Resequencing Barcode {0} (#{1} of {2})'.format(bc, i+1, len(bc_list)))
            output_dir = self.get_output_folder( bc )

            # Extract any consensus sequences associated with this barcode
            record_list = [r for r in amp_analysis_records if r.barcode == bc]
            log.info('Identified {0} consensus sequences for Barcode {1}'.format(len(record_list), bc))
            filtered_records = [r for r in record_list if r.num_reads >= 20]
            unique_records = get_unique_records( filtered_records )
            fraction = 100 * round(len(unique_records)/float(len(record_list)), 3)
            log.info('{0} of {1} ({2}%) consensus sequences passed all filters'.format(len(unique_records),
                                                                                       len(record_list),
                                                                                       fraction))
            reference_file = os.path.join( output_dir, 'reference.fasta' )
            write_records( unique_records, reference_file )

            # Identify all high-quality, barcode-specific ZMWs and write them to file
            zmw_list = get_barcode_zmws( bc_reader, bc )
            zmw_list = filter_zmw_list( bash5, zmw_list, min_snr=min_snr )
            whitelist_file = os.path.join( output_dir, 'whitelist.txt' )
            write_zmw_whitelist( bash5, zmw_list, whitelist_file )

            # Resequence the selected consensus sequences with the selected ZMWs
            self.resequencer( baxh5_file,
                              whitelist_file,
                              reference_file,
                              output=output_dir,
                              min_length=min_length )

            log.info("Finished resequencing Barcode {0}\n".format( bc ))