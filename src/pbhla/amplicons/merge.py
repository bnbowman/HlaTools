#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import logging
import tempfile

from pbcore.io import FastaWriter, FastqReader, FastaRecord, FastqRecord

from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import fasta_size, read_fasta_dict, write_fasta
from pbhla.utils import read_fastq_dict, consensus_filetype, write_fastq
from pbhla.io.BlasrIO import BlasrReader

log = logging.getLogger()
log.setLevel( logging.INFO )

def merge_amplicons( sequence_5p, sequence_3p, output ):
    file_list = [sequence_5p, sequence_3p]
    filetype = consensus_filetype( file_list )
    alignment_file = align_amplicons( filetype, sequence_5p, sequence_3p )
    positions = parse_alignment_positions( alignment_file )
    sequences = read_sequences( filetype, file_list )
    merged = merge_sequences( filetype, sequences, positions )
    write_sequences( filetype, merged, output )

def write_sequences( filetype, merged, output ):
    if filetype == 'fasta':
        write_fasta( merged, output )
    elif filetype == 'fastq':
        write_fastq( merged, output )
    else:
        raise ValueError

def merge_sequences( filetype, sequences, positions ):
    if filetype == 'fastq':
        return merge_fastq_sequences( sequences, positions )
    elif filetype == 'fasta':
        return merge_fasta_sequences( sequences, positions )
    else:
        raise ValueError

def merge_fasta_sequences( sequences, positions ):
    merged = []
    for i, part_list in enumerate( positions ):
        name = "Merged%s_NumReads100" % (i+1)
        sequence = ''
        for part in part_list:
            source = part['name']
            start = part['start']
            end = part['end']
            sequence += sequences[source].sequence[start:end]
        record = FastaRecord(name=name, sequence=sequence)
        merged.append(record)
    return merged

def merge_fastq_sequences( sequences, positions ):
    merged = []
    for i, part_list in enumerate( positions ):
        name = "Merged%s_NumReads100" % (i+1)
        sequence = ''
        quality = ''
        for part in part_list:
            source = part['name']
            start = part['start']
            end = part['end']
            sequence += sequences[source].sequence[start:end]
            quality += sequences[source].qualityString[start:end]
        record = FastqRecord(name=name, sequence=sequence, qualityString=quality)
        merged.append(record)
    return merged

def read_sequences( filetype, file_list ):
    if filetype == 'fastq':
        return read_fastq_dict( file_list )
    elif filetype == 'fasta':
        return read_fasta_dict( file_list )
    else:
        raise ValueError

def parse_alignment_positions( alignment_file ):
    positions = []
    for hit in BlasrReader( alignment_file ):
        left = {'name':hit.qname,
                'start':1,
                'end':int(hit.qstart)}
        right = {'name':hit.tname,
                 'start':int(hit.tstart),
                 'end':int(hit.tlength)}
        positions.append( (left, right) )
    return positions

def write_temp_fasta( fastq_file ):
    """
    Write a temporary Fasta file from a Fastq
    """
    temp = tempfile.NamedTemporaryFile( suffix='.fasta', delete=False )
    with FastaWriter( temp.name ) as handle:
        for record in FastqReader( fastq_file ):
            temp_record = FastaRecord( record.name, record.sequence )
            handle.writeRecord( temp_record )
    return temp

def align_amplicons( filetype, sequence_5p, sequence_3p ):
    blasr_args = {'bestn': 1,
                  'out': 'test.m5',
                  'm': 5,
                  'noSplitSubreads': True}
    if filetype == 'fastq':
        temp_5p = write_temp_fasta( sequence_5p )
        temp_3p = write_temp_fasta( sequence_3p )
        align_left = run_blasr( temp_5p.name, temp_3p.name, blasr_args, verbose=True )
    elif filetype == 'fasta':
        assert fasta_size( sequence_5p ) == 2
        assert fasta_size( sequence_3p ) == 2
        align_left = run_blasr( sequence_5p, sequence_3p, blasr_args )
    else:
        raise ValueError
    return align_left

if __name__ == '__main__':
    import sys

    sequence_5p = sys.argv[1]
    sequence_3p = sys.argv[2]
    output = sys.argv[3]

    logging.basicConfig( stream=sys.stdout )
    merge_amplicons( sequence_5p, sequence_3p, output )