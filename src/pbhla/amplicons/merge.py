#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import is_fasta, write_temp_fasta
from pbhla.utils import is_fastq

def merge_amplicons( sequence_5p, sequence_3p ):
    alignments = align_amplicons( sequence_5p, sequence_3p )
    with open(alignments) as handle:
        for line in handle:
            print line.strip()

def align_amplicons( sequence_5p, sequence_3p ):
    blasr_args = {'bestn': 1}
    if is_fastq( sequence_5p ) and is_fastq( sequence_3p ):
        temp_5p = write_temp_fasta( sequence_5p )
        temp_3p = write_temp_fasta( sequence_3p )
        align_left = run_blasr( temp_5p.name, temp_3p.name, blasr_args )
    elif is_fasta( sequence_5p ) and is_fasta( sequence_3p ):
        align_left = run_blasr( sequence_5p, sequence_3p, blasr_args )
    else:
        raise ValueError
    return align_left

if __name__ == '__main__':
    import sys

    sequence_5p = sys.argv[1]
    sequence_3p = sys.argv[2]

    merge_amplicons( sequence_5p, sequence_3p )