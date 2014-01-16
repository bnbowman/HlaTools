#! /usr/bin/env python

import sys, os, logging, tempfile

from pbhla.fasta.utils import fasta_size, write_temp_fasta
from pbhla.sequences.utils import read_sequences
from pbhla.external.commandline_tools import run_blasr
from pbhla.io.BlasrIO import BlasrReader, BlasrWriter, BlasrM1, BlasrM5, pctsimilarity
from pbhla.utils import check_output_file

NPROC = 4

log = logging.getLogger()

def align_by_identity( query, reference_fasta, output=None, format='1' ):
    """
    Type sequences in a fasta file by finding the closet reference
    """
    # If output isn't specified, base it on the query
    assert format in ['1', '5']
    if output is None:
        basename = '.'.join( query.split('.')[:-1] )
        output = '%s.m%s' % (basename, format)
    # Iterate over each Fasta, aligning individually.
    with BlasrWriter( output ) as handle:
        handle.write_header( 'm1' )
        for record in read_sequences( query ):
            temp = write_temp_fasta( record )
            alignments = _align_fasta( temp.name, reference_fasta, format )
            if not alignments:
                log.info("No hits found for %s" % record.name)
                continue
            alignments = _sort_alignments( alignments )
            alignments = _filter_alignments( alignments )
            #for alignment in alignments:
            #    print alignment
            #    print pctsimilarity( alignment )
            handle.write( alignments[0] )
            os.unlink( temp.name )
    check_output_file( output )
    return output

def _align_fasta( query, reference, format ):
    """
    Align a single query sequence to all valid references
    """
    suffix = '.m%s' % format
    temp_align = tempfile.NamedTemporaryFile( suffix=suffix, delete=False )
    reference_count = fasta_size( reference )
    blasr_args = {'nproc': NPROC,
                  'out': temp_align.name,
                  'bestn': reference_count,
                  'nCandidates': reference_count,
                  'm': format,
                  'noSplitSubreads': True}
    run_blasr( query, reference, blasr_args )
    # Parse the output for return and delete the file
    alignments = list( BlasrReader( temp_align.name ))
    os.unlink( temp_align.name )
    return alignments

def _sort_alignments( alignments ):
    """
    Sort alignments by Percent-Identity first, and Score second
    """
    alignments = sorted( alignments, key=lambda x: pctsimilarity(x), reverse=True )
    return alignments

def _filter_alignments( alignments ):
    """
    Filter out all but the best matches
    """
    max_pctid = pctsimilarity(alignments[0])
    alignments = filter( lambda x: pctsimilarity(x) >= max_pctid,
                         alignments )
    return alignments

if __name__ == '__main__':
    logging.basicConfig( level=logging.INFO )

    query_fasta = sys.argv[1]
    reference_fasta = sys.argv[2]
    output_file = sys.argv[3]
    format = sys.argv[4]

    align_by_identity( query_fasta, reference_fasta, output_file, format )
