#! /usr/bin/env python

import os, logging, tempfile

from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter
from pbhla.external.commandline_tools import run_blasr
from pbhla.utils import check_output_file, remove_file, read_list_file
from pbhla.fasta.utils import write_fasta
from pbhla.io.BlasrIO import BlasrReader

NPROC = 4
TEMP_FASTA = 'temp.fasta' 

log = logging.getLogger()

def extract_exons( input_file, exon_fofn ):
    """
    Extract all exons from a particular Fasta File into a separate Fasta File
    """
    basename = '.'.join( input_file.split('.')[:-1] )
    output_file = '%s.exons.fasta' % basename
    fasta_record = _read_fasta_record( input_file )
    with FastaWriter( output_file ) as output:
        for exon in _extract_exons( fasta_record, exon_fofn, basename ):
            output.writeRecord( exon )
    return output_file

def _read_fasta_record( input_file ):
    records = list(FastaReader( input_file ))
    if len( records ) == 1:
        return records[0]
    msg = 'expected a single Fasta, found MultiFasta!'
    log.error( msg )
    raise TypeError( msg )

def _extract_exons( fasta_record, exon_fofn, prefix ):
    """
    Extract Fasta Records of each Exon from a Fasta Record
    """
    for exon_fasta in read_list_file( exon_fofn ):
        exon_num = exon_fasta[-7]
        start, end = _find_exon_position( fasta_record, exon_fasta, prefix )
        if start is None or end is None:
            log.info('No Exon #%s to extract from "%s"' % (exon_num, fasta_record.name))
            continue
        yield _extract_exon_sequence( fasta_record, exon_num, start, end )

def _find_exon_position( fasta_record, exon_fasta, prefix ):
    """
    Find the start and end position of a certain Exon in a given Fasta Record
    """
    temp = tempfile.NamedTemporaryFile( suffix='.fasta', delete=False )
    write_fasta( [fasta_record], temp.name )
    alignment_file = _align_exons( temp.name, exon_fasta, prefix )
    os.unlink( temp.name )
    if alignment_file:
        return _parse_exon_location( alignment_file )
    return None, None

def _align_exons( query_fasta, exon_fasta, prefix ):
    """
    Align all supplied exon sequences to a "Query" Fasta
    """
    exon_num = exon_fasta[-7]
    alignment_file = '%s.exon%s.m1' % (prefix, exon_num)
    log.info('Attempting')
    blasr_args = {'nproc': NPROC,
                  'out': alignment_file,
                  'bestn': 1,
                  'noSplitSubreads': True}
    run_blasr( exon_fasta,
               query_fasta,
               blasr_args )
    count = count_hits( alignment_file )
    if count > 0:
        log.info('Found %s hits for Exon #%s' % (count, exon_num))
        return alignment_file
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

def _extract_exon_sequence( fasta, exon_num, start, end ):
    """
    Create an Exon record from its coordinates and a Fasta
    """
    exon_name = '%s_exon%s' % (fasta.name, exon_num)
    exon_sequence = fasta.sequence[start:end]
    return FastaRecord( exon_name, exon_sequence )

### Utilities ###

def count_hits( filename ):
    count = 0
    with open(filename, 'r') as handle:
        for line in handle:
            count += 1
    return count

if __name__ == '__main__':
    import sys
    logging.basicConfig( stream=sys.stdout, level=logging.INFO )

    input_fasta = sys.argv[1]
    exon_fofn = sys.argv[2]

    extract_all_exons( input_fasta, exon_fofn )
