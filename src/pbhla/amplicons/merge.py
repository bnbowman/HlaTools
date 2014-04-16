#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import logging

from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import is_fasta, write_temp_fasta

log = logging.getLogger()

def is_fastq( filename ):
    if filename.endswith('.fastq') or filename.endswith('.fq'):
        return True
    return False

def merge_amplicons( sequence_5p, sequence_3p ):
    alignments = align_amplicons( sequence_5p, sequence_3p )
    #with open(alignments) as handle:
    #    for line in handle:
    #        print line.strip()

def align_amplicons( sequence_5p, sequence_3p ):
    blasr_args = {'bestn': 1, 'out': 'test.m5', 'm': 5}
    if is_fastq( sequence_5p ) and is_fastq( sequence_3p ):
        temp_5p = write_temp_fasta( sequence_5p )
        temp_3p = write_temp_fasta( sequence_3p )
        print temp_5p, temp_5p.name
        print temp_3p, temp_3p.name
        align_left = run_blasr( temp_5p.name, temp_3p.name, blasr_args, verbose=True )
    elif is_fasta( sequence_5p ) and is_fasta( sequence_3p ):
        print sequence_5p
        print sequence_3p
        align_left = run_blasr( sequence_5p, sequence_3p, blasr_args )
    else:
        raise ValueError
    print align_left
    return align_left

if __name__ == '__main__':
    import sys

    sequence_5p = sys.argv[1]
    sequence_3p = sys.argv[2]

    logging.basicConfig( stream=sys.stdout, level=logging.INFO )
    merge_amplicons( sequence_5p, sequence_3p )