#! /usr/bin/env python
import os
from pbhla.utils import log

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'


def get_file_type( filename ):
    if filename.endswith('.fa') or filename.endswith('.fasta'):
        return 'fasta'
    elif filename.endswith('.fq') or filename.endswith('.fastq'):
        return 'fastq'
    elif filename.endswith('.fofn'):
        return 'fofn'
    elif filename.endswith('.bas.h5') or filename.endswith('.bax.h5'):
        return 'bas.h5'
    else:
        msg = 'File is not of a recognized filetype'
        log.error( msg )
        raise TypeError( msg )

def get_file_source( filename ):
    base_name = os.path.basename( filename )
    root_name = base_name.split('.')[0]
    parts = root_name.split('_')
    return parts[1]

def is_fasta( filename ):
    if filename.endswith('.fasta') or filename.endswith('.fa'):
        return True
    return False

def is_fastq( filename ):
    if filename.endswith('.fastq') or filename.endswith('.fq'):
        return True
    return False

def get_file_root( filename ):
    return '.'.join( filename.split('.')[:-1] )