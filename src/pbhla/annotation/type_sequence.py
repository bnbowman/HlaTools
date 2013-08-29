#! /usr/bin/env python

import sys, os, logging

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbhla.fasta.utils import write_fasta, fasta_size
from pbhla.external.commandline_tools import run_blasr
from pbhla.io.BlasrIO import BlasrReader, BlasrWriter

TMP_FASTA = 'temp.fasta'
TMP_ALIGN = 'temp.m1'
NPROC = 4

log = logging.getLogger()

def align_by_identity( query_fasta, reference_fasta, output_file ):
    """
    Type sequences in a fasta file by finding the closet reference
    """
    with BlasrWriter( output_file ) as handle:
        handle.writeHeader()
        for record in FastaReader( query_fasta ):
            write_fasta( record, TMP_FASTA )
            alignments = _align_fasta( TMP_FASTA, reference_fasta )
            alignments = _sort_alignments( alignments )
            alignments = _filter_alignments( alignments )
            for alignment in alignments:
                handle.write( alignment )
            os.remove( TMP_FASTA )

def _align_fasta( query, reference ):
    """
    Align a single query sequence to all valid references
    """
    reference_count = fasta_size( reference )
    blasr_args = {'nproc': NPROC,
                  'out': TMP_ALIGN,
                  'bestn': reference_count,
                  'nCandidates': reference_count,
                  'noSplitSubreads': True}
    run_blasr( query, reference, blasr_args )
    # Parse the output for return and delete the file
    alignments = list( BlasrReader( TMP_ALIGN ))
    os.remove( TMP_ALIGN )
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
