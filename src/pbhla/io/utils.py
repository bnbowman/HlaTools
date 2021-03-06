#! /usr/bin/env python
from pbhla.filenames import get_file_type

__author__ = 'bbowman@pacificbiosciences.com'

import logging
from tempfile import NamedTemporaryFile

from pbcore.io.FastaIO import FastaRecord, FastaWriter
from pbcore.io.FastqIO import FastqRecord, FastqWriter
from pbhla.io.AmpAnalysisIO import AmpliconAnalysisRecord, AmpliconAnalysisWriter
from pbhla.utils import check_output_file

log = logging.getLogger()

def slice_2d( args ):
    """
    Convert a __getitems__ input into a pair of slice objects for ease-of-use
    """
    if isinstance( args, int ):
        return _to_slice(args), slice(None)
    elif isinstance( args, slice ):
        return args, slice(None)
    elif isinstance( args, tuple ):
        if len(args) > 2:
            raise ValueError("Cannot create 2D slice from %s arguments" % len(args))
        else:
            return (_to_slice(item) for item in args)
    else:
        raise TypeError("slice_2d accepts Int, Slice, or Tuple arguments")

def _to_slice( item ):
    """
    Convert an item from Int to Slice if needed
    """
    if isinstance( item, int ):
        return slice( item, item+1 )
    if isinstance( item, slice ):
        return item
    else:
        raise TypeError("Input must be Int or Slice")

def parse_locus_dict( filename ):
    """
    Read a dictionary of values from a file with locus-specific data
    """
    data = {}
    with open( filename, 'r' ) as handle:
        for line in handle:
            datum, locus = line.strip().split()
            if locus in data:
                msg = 'Duplicate locus fofn "%s"!' % locus
                log.error( msg )
                raise KeyError( msg )
            else:
                data[locus] = datum
    return data

def get_output_file( input_file, modifier ):
    """
    Get a modified output file name based on some input file
    """
    basename = '.'.join( input_file.split('.')[:-1] )
    file_type = get_file_type( input_file )
    return '%s.%s.%s' % (basename, modifier, file_type)

def get_temp_fasta_record( record ):
    """
    If a record isn't in Fasta format, try to create a FastaRecord from it
    """
    if isinstance( record, FastaRecord ):
        return record
    try:
        return FastaRecord( record.name.strip(), record.sequence.strip() )
    except:
        msg = 'Unrecognized sequence record type'
        log.error( msg )
        raise TypeError( msg )

def get_temp_fasta( record ):
    """
    Create a temporary Fasta file for Blasr/HMMsearch/etc
    """
    temp_record = get_temp_fasta_record( record )
    temp_fasta = NamedTemporaryFile( suffix='.fasta' )
    with FastaWriter( temp_fasta.name ) as handle:
        handle.writeRecord( temp_record )
    return temp_fasta

def write_records( records, filename ):
    if all([isinstance(r, FastaRecord) for r in records]):
        write_fasta_records( records, filename )
    elif all([isinstance(r, FastqRecord) for r in records]):
        write_fastq_records( records, filename )
    elif all([isinstance(r, AmpliconAnalysisRecord) for r in records]):
        write_amp_analysis_records( records, filename )
    else:
        msg = 'Invalid sequence record type'
        log.error( msg )
        raise TypeError( msg )

def write_fasta_records( records, filename ):
    log.info("Writing {0} FastaRecords to {1}".format(len(records), filename))
    with FastaWriter( filename ) as handle:
        for record in records:
            handle.writeRecord( record )
    check_output_file( filename )

def write_fastq_records( records, filename ):
    log.info("Writing {0} FastqRecords to {1}".format(len(records), filename))
    with FastqWriter( filename ) as handle:
        for record in records:
            handle.writeRecord( record )
    check_output_file( filename )

def write_amp_analysis_records( records, filename ):
    log.info("Writing {0} AmpAnalysisRecords to {1}".format(len(records), filename))
    with AmpliconAnalysisWriter( filename ) as handle:
        for record in records:
            handle.write_fasta( record )
    check_output_file( filename )

def get_unique_records( records ):
    """Remove redundant sequences, primarily to avoid confusing Quiver"""
    sorted_records = sorted(records, key=lambda r: len(r.sequence), reverse=True)

    unique_sequences = []
    unique_records = []
    for record in sorted_records:
        is_unique = True
        for sequence in unique_sequences:
            if record.sequence in sequence:
                is_unique = False
                break
        if is_unique:
            unique_records.append( record )
            unique_sequences.append( record.sequence )

    return unique_records
