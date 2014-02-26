#! /usr/bin/env python

import os 
import logging
from operator import itemgetter

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbcore.io.FastqIO import FastqReader, FastqWriter
from pbhla.fasta.utils import fasta_size, write_fasta
from pbhla.external.utils import get_alignment_file
from pbhla.filenames import get_file_type
from pbhla.io.BlasrIO import BlasrReader
from pbhla.utils import check_output_file

log = logging.getLogger()

def extract_alleles( input_file, min_reads, min_length, output_file=None ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    # Set the output file if not specified
    output_file = output_file or _get_output_file( input_file )
    output_type = get_file_type( output_file )
    # Parse the alignment data and extract the target sequences
    sequences = _parse_input_records( input_file )
    sequences = _filter_on_length( sequences, min_length )
    sequences = _filter_on_numreads( sequences, min_reads )
    _write_output( sequences, output_file, output_type )
    return output_file

def _get_output_file( input_file ):
    basename = '.'.join( input_file.split('.')[:-1] )
    file_type = get_file_type( input_file )
    return '%s.filtered.%s' % (basename, file_type)

def _parse_input_records( input_file ):
    """
    Parse the input sequence records with the appropriate pbcore Reader
    """
    input_type = get_file_type( input_file )
    if input_type == 'fasta':
        return list( FastaReader( input_file ))
    elif input_type == 'fastq':
        return list( FastqReader( input_file ))
    else:
        msg = 'Input file must be either Fasta or Fastq'
        log.error( msg )
        raise TypeError( msg )

def _filter_on_numreads( sequences, min_reads ):
    """
    Filter reads based on the NumReads from AmpAnalysis
    """
    # Internal function for simplicity
    def record_size( record ):
        return int( record.name.split('NumReads')[-1] )
    return [r for r in sequences if record_size( r ) >= min_reads]

def _filter_on_length( sequences, min_length ):
    """
    Filter out reads below a certain length
    """
    return [r for r in sequences if len(r.sequence) >= min_length]

def _write_output( records, output_file, output_type ):
    """
    Write the records out to file
    """
    if output_type == 'fasta':
        write_fasta( records, output_file )
    else:
        with FastqWriter( output_file ) as writer:
            for record in records:
                writer.writeRecord( record )
    check_output_file( output_file )
    return output_file

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_file = sys.argv[1]
    min_reads = int(sys.argv[2])
    min_length = int(sys.argv[3])
    
    extract_alleles( input_file, min_reads, min_length )
