#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

from pbcore.io.BarcodeH5Reader import BarcodeH5Fofn, BarcodeH5Reader
from pbhla.utils import validate_file

def get_barcode_reader( barcode_file ):
    barcode_file = validate_file( barcode_file )
    if barcode_file.endswith('.fofn'):
        return BarcodeH5Fofn( barcode_file )
    elif barcode_file.endswith('.bc.h5'):
        return BarcodeH5Reader( barcode_file )
    else:
        raise ValueError

def barcode_string_to_list( barcode_string ):
    if barcode_string is None:
        return None
    barcodes = barcode_string.split(',')
    return [format_barcode(b) for b in barcodes]

def format_barcode( barcode ):
    if '--' in barcode:
        return barcode
    else:
        return 'F{0}--R{0}'.format( barcode )

def get_barcodes( bc_reader, bc_string ):
    barcodes = barcode_string_to_list( bc_string )
    for label in bc_reader.barcodeLabels:
        print barcodes, label