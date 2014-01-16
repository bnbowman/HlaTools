#! /usr/bin/env python

import os, logging, tempfile

from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter
from pbcore.io.FastqIO import FastqRecord, FastqWriter

from pbhla.external.commandline_tools import run_blasr
from pbhla.utils import (check_output_file, 
                         remove_file, 
                         read_list_file,
                         get_file_type,
                         count_hits)
from pbhla.fasta.utils import write_fasta
from pbhla.io.BlasrIO import BlasrReader

NPROC = 4
TEMP_FASTA = 'temp.fasta' 

log = logging.getLogger()

def extract_exons( input_record, exon_fofn, directory=None ):
    """
    Extract all exons from a particular Fasta File into a separate Fasta File
    """
    if isinstance( input_record, str ):
        output_type = get_file_type( input_record )
        input_record = _read_fasta_record( input_record )
    elif isinstance( input_record, FastaRecord ):
        output_type = 'fasta'
    elif isinstance( input_record, FastqRecord ):
        output_type = 'fastq'
    else:
        msg = 'Input record must be Filename, FastaRecord or FastqRecord'
        log.error( msg )
        raise TypeError( msg )
    return _extract_exons( input_record, exon_fofn, output_type, directory )

def _read_fasta_record( input_file ):
    """
    Read a single FastaRecord, raising an error if a MultiFasta is found
    """
    records = list(FastaReader( input_file ))
    if len( records ) == 1:
        return records[0]
    msg = 'expected a single Fasta, found MultiFasta!'
    log.error( msg )
    raise TypeError( msg )

def _extract_exons( record, exon_fofn, output_type, directory ):
    """
    Extract Fasta Records of each Exon from a Fasta Record
    """
    log.info('Extracting exons from "%s"' % record.name)
    temp_fasta = _write_temp_fasta( record )
    output_file = os.path.join( directory, 'all_exons.%s' % output_type )
    output_handle = _open_output_handle( output_file, output_type )

    # Iterate over the individual Exon Fasta files looking for alignments
    exon_count = 0
    for exon_fasta in read_list_file( exon_fofn ):
        exon_num = exon_fasta[-7]
        start, end = _find_exon_position( temp_fasta, exon_fasta, directory )
        if start is None or end is None:
            continue
        exon_count += 1
        exon_record = _extract_exon_record( record, exon_num, start, end )
        output_handle.writeRecord( exon_record )
    os.unlink( temp_fasta )
    output_handle.close()

    if exon_count:
        log.info("Extracted %s exons from %s" % (exon_count, record.name))
    else:
        log.warn("No valid exons found for %s!" % record.name)
        return None
    check_output_file( output_file )
    return output_file

def _open_output_handle( output_file, output_type ):
    """
    Open an appropriate output handle to record the exon sequences
    """
    if output_type == 'fasta':
        return FastaWriter( output_file )
    elif output_type == 'fastq':
        return FastqWriter( output_file )
    msg = 'Output type must be Fasta or Fastq'
    log.error( msg )
    raise TypeError( msg )

def _write_temp_fasta( record ):
    """
    Write a sequence record out to a temporary Fasta file
    """
    temp = tempfile.NamedTemporaryFile( suffix='.fasta', delete=False )
    if isinstance( record, FastaRecord ):
        write_fasta( [record], temp.name )
    elif isinstance( record, FastqRecord ):
        temp_record = FastaRecord(record.name, record.sequence)
        write_fasta( [temp_record], temp.name )
    else:
        msg = 'Record must be either FastaRecord or FastqRecord'
        log.error( msg )
        raise TypeError( msg )
    return temp.name

def _find_exon_position( query, exon_fasta, directory ):
    """
    Find the start and end position of a certain Exon in a given Fasta Record
    """
    alignment_file = _align_exons( query, exon_fasta, directory )
    if alignment_file:
        return _parse_exon_location( alignment_file )
    return None, None

def _align_exons( query, exon_fasta, directory ):
    """
    Align all supplied exon sequences to a "Query" Fasta
    """
    exon_num = exon_fasta[-7]
    alignment_file = 'exon%s.m1' % exon_num
    alignment_path = os.path.join( directory, alignment_file )
    blasr_args = {'nproc': NPROC,
                  'out': alignment_path,
                  'minSubreadLength': 30,
                  'minReadLength': 30,
                  'maxScore': 0,
                  'bestn': 1,
                  'noSplitSubreads': True}
    run_blasr( exon_fasta,
               query,
               blasr_args )
    count = count_hits( alignment_path )
    if count > 0:
        log.info('Found %s hits for Exon #%s' % (count, exon_num))
        return alignment_path
    log.info('No hits found for Exon #%s' % exon_num)
    return None

def _parse_exon_location( alignment_file ):
    """
    Parse the most likely Exon location from an Exon-Fasta alignment
    """
    alignments = list( BlasrReader( alignment_file ))
    alignments = sorted( alignments, key=lambda x: int(x.score))
    alignments = sorted( alignments, key=lambda x: float(x.pctsimilarity), reverse=True)
    return int(alignments[0].tstart), int(alignments[0].tend)

def _extract_exon_record( record, exon_num, start, end ):
    """
    Create an Exon record from its coordinates and a Fasta
    """
    exon_name = '%s_exon%s' % (record.name, exon_num)
    exon_sequence = record.sequence[start:end]
    if isinstance( record, FastaRecord ):
        return FastaRecord( exon_name, exon_sequence )
    elif isinstance( record, FastqRecord ):
        exon_qual = record.qualityString[start:end]
        return FastqRecord( exon_name, exon_sequence, qualityString=exon_qual)
    msg = 'Record must be either FastaRecord or FastqRecord'
    log.error( msg )
    raise TypeError( msg )

if __name__ == '__main__':
    import sys
    logging.basicConfig( stream=sys.stdout, level=logging.INFO )

    input_fasta = sys.argv[1]
    exon_fofn = sys.argv[2]

    extract_all_exons( input_fasta, exon_fofn )
