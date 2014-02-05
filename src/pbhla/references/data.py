#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
import logging, logging.config

from pbcore.io import FastaReader, FastaWriter

from pbhla import __LOG__, __DATA__
from pbhla.utils import valid_file

logging.config.fileConfig( __LOG__ )
log = logging.getLogger(__name__)

GEN_PATTERN  = '_gen.fasta'
CDNA_PATTERN = '_nuc.fasta'
GEN_REF  = os.path.join( __DATA__, 'genomic.fasta')
CDNA_REF = os.path.join( __DATA__, 'cDNA.fasta')

def get_genomic_reference():
    if valid_file( GEN_REF ):
        return GEN_REF
    else:
        create_genomic_reference()
        return get_genomic_reference()

def get_cDNA_reference():
    if valid_file( CDNA_REF ):
        return CDNA_REF
    else:
        create_cDNA_reference()
        return get_cDNA_reference()

def list_data_files( end_pattern ):
    for root, dirs, files in os.walk( __DATA__ ):
        for filename in files:
            if filename.endswith( end_pattern ):
                yield os.path.join(root, filename)

def append_records( handle, filepath ):
    for record in FastaReader( filepath ):
        handle.writeRecord( record )

def create_genomic_reference():
    with FastaWriter( GEN_REF ) as handle:
        for filepath in list_data_files( GEN_PATTERN ):
            append_records( handle, filepath )
    return GEN_REF

def create_cDNA_reference():
    with FastaWriter( CDNA_REF ) as handle:
        for filepath in list_data_files( CDNA_PATTERN ):
            append_records( handle, filepath )
    return CDNA_REF

if __name__ == '__main__':
    print get_genomic_reference()
    print get_cDNA_reference()