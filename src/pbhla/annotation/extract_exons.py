#! /usr/bin/env python

import sys, os, logging

from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter
from pbhla.external.commandline_tools import run_blasr
from pbhla.utils import check_output_file
from pbhla.fasta.utils import write_fasta
from pbhla.io.BlasrIO import BlasrReader

NPROC = 4
TEMP_ALIGN = 'temp.m1'
TEMP_FASTA = 'temp.fasta' 

log = logging.getLogger()

def extract_all_exons( input_fasta, exon_fofn ):
    """
    Extract all exons from a particular Fasta File into a separate Fasta File
    """
    fasta_record = list(FastaReader( input_fasta ))[0]
    output_file = _get_output_file( input_fasta )
    with FastaWriter( output_file ) as output:
        for exon in _extract_exons( fasta_record, exon_fofn ):
            output.writeRecord( exon )

def _extract_exons( fasta_record, exon_fofn ):
    """
    Extract Fasta Records of each Exon from a Fasta Record
    """
    for exon_fasta in _parse_fofn( exon_fofn ):
        exon_num = exon_fasta[-7]
        start, end = _find_exon_position( fasta_record, exon_fasta )
        if start is None or end is None:
            log.info('No Exon #%s to extract from "%s"' % (exon_num, fasta_record.name))
            continue
        yield _extract_exon_sequence( fasta_record, exon_num, start, end )

def _find_exon_position( fasta_record, exon_fasta ):
    """
    Find the start and end position of a certain Exon in a given Fasta Record
    """
    write_fasta( [fasta_record], TEMP_FASTA )
    alignment_file = _align_exons( TEMP_FASTA, exon_fasta )
    remove_file( TEMP_FASTA )
    if alignment_file:
        return _parse_exon_location( alignment_file )
    return None, None

def _align_exons( query_fasta, exon_fasta ):
    """
    Align all supplied exon sequences to a "Query" Fasta
    """
    remove_file( TEMP_ALIGN )
    exon_num = exon_fasta[-7]
    log.info('Attempting')
    blasr_args = {'nproc': NPROC,
                  'out': TEMP_ALIGN,
                  'bestn': 1,
                  'noSplitSubreads': True}
    run_blasr( exon_fasta,
               query_fasta,
               blasr_args )
    count = count_hits( TEMP_ALIGN )
    if count > 0:
        log.info('Found %s hits for Exon #%s' % (count, exon_num))
        return TEMP_ALIGN
    log.info('No hits found for Exon #%s' % exon_num)
    return None

def _parse_exon_location( alignment_file ):
    """
    Parse the most likely Exon location from an Exon-Fasta alignment
    """
    alignments = list( BlasrReader( alignment_file ))
    alignments = sorted( alignments, key=lambda x: int(x.score))
    alignments = sorted( alignments, key=lambda x: float(x.pctsimilarity), reverse=True)
    print alignments
    print alignments[0]
    print alignments[-1]
    return int(alignments[0].tstart), int(alignments[0].tend)

def _extract_exon_sequence( fasta, exon_num, start, end ):
    """
    Create an Exon record from its coordinates and a Fasta
    """
    exon_name = '%s_exon%s' % (fasta.name, exon_num)
    exon_sequence = fasta.sequence[start:end]
    return FastaRecord( exon_name, exon_sequence )

### Utilities ###

def _get_output_file( filename ):
    basename = '.'.join( filename.split('.')[:-1] )
    return '%s_exons.fasta' % basename

def _parse_fofn( exon_fofn ):
    with open(exon_fofn, 'r') as handle:
        for line in handle:
            yield line.strip()

def count_hits( filename ):
    count = 0
    with open(filename, 'r') as handle:
        for line in handle:
            count += 1
    return count

def remove_file( filename ):
    try:
        os.remove( filename )
    except:
        pass

if __name__ == '__main__':
    input_fasta = sys.argv[1]
    exon_fofn = sys.argv[2]

    logging.basicConfig( stream=sys.stdout, level=logging.INFO )

    extract_all_exons( input_fasta, exon_fofn )
