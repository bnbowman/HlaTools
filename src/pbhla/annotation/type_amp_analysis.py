#! /usr/bin/env python

import os

from pbcore.io import FastaReader, FastaWriter


from pbhla.fasta.orient_fasta import orient_fasta

from pbhla.annotation.top_amp_analysis import top_amp_analysis

NPROC = 6

def type_amp_analysis( input_file, reference_file, output_file ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    alignments = _align_sequences( input_file, reference_file )
    groups = _group_by_locus( alignments )
    ordered = _sort_groups( groups )
    selected = list( _select_sequences( ordered ))
    sequences = list( FastaReader( input_file ))
    subset = list( _subset_sequences( sequences, selected ))
    _write_sequences( subset, output_file )

def _align_sequences( query, reference ):
    """
    Align the output of AA to the references and return
    """
    dirname = os.path.dirname( query )
    alignment_file = os.path.join( dirname, 'temp.m1' )
    # Run Blasr
    ref_count = fasta_size( reference )
    blasr_args = {'nproc': NPROC,
                  'out': alignment_file,
                  'bestn': 1,
                  'nCandidates': ref_count,
                  'noSplitSubreads': True }
    run_blasr( query, reference, blasr_args )
    # Parse and return the Blasr records
    alignments = list( BlasrReader( alignment_file ))
    os.remove( alignment_file )
    return alignments

if __name__ == '__main__':
    import sys

    input_file = sys.argv[1]
    reference_file = sys.argv[2]
    output_file = sys.argv[3]
    
    top_amp_analysis( input_file, reference_file, output_file )
