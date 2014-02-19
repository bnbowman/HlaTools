#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os

from pbhla.resequencing.Resequencer import Resequencer
from pbhla.barcode.utils import get_barcode_reader, get_barcodes
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

    def __call__(self, data_file, barcode_file, amp_analysis, output=None, barcode_string=None):
        data_file = validate_file( data_file )
        bc_reader = get_barcode_reader( barcode_file )

        bc_list = get_barcodes( bc_reader, barcode_string )
        print data_file
        print barcode_file
        print amp_analysis
        print bc_list