#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
import logging, logging.config

from pbhla import __LOG__
from pbhla.resequencing.Resequencer import Resequencer
from pbhla.bash5 import get_bash5_reader, filter_zmw_list
from pbhla.barcodes import get_barcode_reader, get_barcodes, get_barcode_zmws
from pbhla.sequences.input import get_input_file

logging.config.fileConfig( __LOG__ )
log = logging.getLogger( __name__ )

class AmpliconAnalysisResequencer( object ):
    """
    An object for iterating over all of the barcodes in
    """

    def __init__(self, setup=None, nproc=1):
        """Initialize cross-cluster and object-specific settings"""
        self._resequencer = Resequencer(setup, nproc)

    @property
    def resequencer(self):
        return self._resequencer

    @property
    def setup(self):
        return self.resequencer.setup

    @property
    def nproc(self):
        return self.resequencer.nproc

    def __call__(self, data_file, barcode_file, amp_analysis, output=None, barcode_string=None):
        bash5 = get_bash5_reader( data_file )
        bc_reader = get_barcode_reader( barcode_file )
        amp_analysis = get_input_file( amp_analysis )

        bc_list = get_barcodes( bc_reader, barcode_string )
        for i, bc in enumerate( bc_list ):
            log.info('Analyzing Barcode {0} (#{1} of {2})'.format(bc, i+1, len(bc_list)))
            zmw_list = get_barcode_zmws( bc_reader, bc )
            zmw_list = filter_zmw_list( bash5, zmw_list, min_snr=3.0 )