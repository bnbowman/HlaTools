#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import logging

from pbcore.io.BarcodeH5Reader import BarcodeH5Fofn, BarcodeH5Reader
from pbhla.utils import validate_file

log = logging.getLogger()

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

    # If no barcodes specified, return them all
    if barcodes is None:
        return bc_reader.barcodeLabels

    # Otherwise return the intersection and warn if missing
    bc_list = []
    for bc in barcodes:
        if bc in bc_reader.barcodeLabels:
            bc_list.append( bc )
        else:
            log.warn('Barcode "%s" requested but not found, skipping' % bc)
    return bc_list

def get_barcode_reads( bc_reader, bc ):
    if isinstance( bc_reader, BarcodeH5Reader ):
        reads = bc_reader.labeledZmwsFromBarcodeLabel( bc )
    elif isinstance( bc_reader, BarcodeH5Fofn ):
        reads = reduce(lambda x,y: x + y,
                       map(lambda z: z.labeledZmwsFromBarcodeLabel( bc ),
                           bc_reader._bcH5s))
    else:
        raise ValueError
    log.info("Identified {0} reads from Barcode {1}".format(len(reads), bc))
    return reads