#! /usr/bin/env python

import sys, os, logging, tempfile

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbhla.fasta.utils import write_fasta, fasta_size
from pbhla.external.commandline_tools import run_blasr
from pbhla.io.BlasrIO import BlasrReader, BlasrWriter
from pbhla.utils import check_output_file

NPROC = 4

log = logging.getLogger()

def align_by_identity( query_fasta, reference_fasta, output=None ):
    """
    Type sequences in a fasta file by finding the closet reference
    """
    # If output isn't specified, base it on the query
    if output is None:
        basename = '.'.join( query_fasta.split('.')[:-1] )
        output = '%s.m1' % basename
    # Iterate over each Fasta, aligning individually.
    with BlasrWriter( output ) as handle:
        handle.writeHeader()
        for record in FastaReader( query_fasta ):
            temp = tempfile.NamedTemporaryFile( suffix='.fasta', delete=False )
            write_fasta( record, temp.name )
            alignments = _align_fasta( temp.name, reference_fasta )
            alignments = _sort_alignments( alignments )
            alignments = _filter_alignments( alignments )
            for alignment in alignments:
                handle.write( alignment )
            os.unlink( temp.name )
    check_output_file( output )
    return output

def _align_fasta( query, reference ):
    """
    Align a single query sequence to all valid references
    """
    temp_align = tempfile.NamedTemporaryFile( suffix='.m1', delete=False )
    reference_count = fasta_size( reference )
    blasr_args = {'nproc': NPROC,
                  'out': temp_align.name,
                  'bestn': reference_count,
                  'nCandidates': reference_count,
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
    alignments = sorted( alignments, key=lambda x: int(x.score))
    alignments = sorted( alignments, key=lambda x: float(x.pctsimilarity),
                                                   reverse=True)
    return alignments

def _filter_alignments( alignments ):
    """
    Filter out all but the best matches
    """
    max_pctid = float(alignments[0].pctsimilarity)
    alignments = filter( lambda x: float(x.pctsimilarity) >= max_pctid,
                         alignments )
    return alignments

if __name__ == '__main__':
    logging.basicConfig( level=logging.INFO )

    query_fasta = sys.argv[1]
    reference_fasta = sys.argv[2]
    output_file = sys.argv[3]

    align_by_identity( query_fasta, reference_fasta, output_file )
