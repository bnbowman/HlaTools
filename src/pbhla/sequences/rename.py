#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import logging

from pbcore.io import FastaRecord, FastqRecord
from pbhla.sequences.utils import read_sequences, write_sequences

log = logging.getLogger()

def rename_sequences( sequence_file ):
    records = read_sequences( sequence_file )
    if any([r.name.strip().endswith('|quiver') for r in records]):
        records = [rename_record(r) for r in records]
        output_file = get_output_file( sequence_file )
        write_sequences( records, output_file )
        return output_file
    else:
        return sequence_file

def rename_record( record ):
    new_name = '|'.join( record.name.strip().split('|')[:-1] )
    if isinstance( record, FastaRecord ):
        return FastaRecord( new_name, record.sequence )
    elif isinstance( record, FastqRecord ):
        return FastqRecord( new_name, record.sequence, record.quality )
    else:
        msg = "Object must be a valid Fasta or Fastq Record"
        log.error( msg )
        raise TypeError( msg )

def get_output_file( input_file ):
    parts = input_file.split('.')
    root = '.'.join( parts[:-1] )
    suffix = parts[-1]
    return '%s.renamed.%s' % (root, suffix)