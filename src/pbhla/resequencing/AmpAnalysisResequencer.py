#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os

from pbhla.resequencing.Resequencer import Resequencer
from pbhla.utils import validate_file

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

    def __call__(self, data_file, barcode_file, amp_analysis, output=None, barcodes=None):
        data_file = validate_file( data_file )
        barcode_file = validate_file( barcode_file )
        barcodes = get_barcodes( barcodes )
        print data_file
        print barcode_file
        print barcodes

def get_barcodes( barcode_string ):
    if barcode_string is None:
        return None
    barcodes = barcode_string.split(',')
    return [format_barcode(b) for b in barcodes]

def format_barcode( barcode ):
    if '--' in barcode:
        return barcode
    else:
        return 'F{0}--R{0}'.format( barcode )