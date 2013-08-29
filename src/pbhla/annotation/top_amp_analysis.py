#! /usr/bin/env python

import os
from operator import itemgetter

from pbcore.io import FastaReader, FastaWriter
from pbhla.fasta.utils import fasta_size
from pbhla.external.commandline_tools import run_blasr
from pbhla.io.BlasrIO import BlasrReader

NPROC = 6

def top_amp_analysis( input_file, reference_file, output_file ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    alignment_file = align_best_reference( input_file, reference_file )
    alignments = list( BlasrReader( alignment_file ))
    groups = _group_by_locus( alignments )
    ordered = _sort_groups( groups )
    selected = list( _select_sequences( ordered ))
    sequences = list( FastaReader( input_file ))
    subset = list( _subset_sequences( sequences, selected ))
    _write_sequences( subset, output_file )

def _group_by_locus( alignments ):
    """
    Group reads by the locus of their best alignment
    """
    loci = {}
    for record in alignments:
        reference = record.tname.split('*')[0]
        locus = reference.split('_')[-1]
        try:
            loci[locus].append( record )
        except:
            loci[locus] = [record]
    return loci

def _sort_groups( groups ):
    """
    Order each group of records individually
    """
    ordered = {}
    for locus, group in groups.iteritems():
        ordered[locus] = _sort_group( group )
    return ordered

def _sort_group( group ):
    """
    Order records in a group by their number of reads
    """
    # Internal function for simplicity
    def record_size( record ):
        return int( record.qname.split('NumReads')[-1] )
    # Count, sort and return
    counts = {record: record_size(record) for record in group}
    tuples = sorted( counts.iteritems(), key=itemgetter(1), reverse=True )
    return [t[0] for t in tuples]

def _select_sequences( groups ):
    """
    Select the top 1-2 sequences for each Locus
    """
    for group in groups.itervalues():
        first, the_rest = group[0], group[1:]
        # Yeild the first sequence from each group
        yield first.qname
        # Yeild the second, if any pass filter
        for record in the_rest:
            if record.tname != first.tname:
                yield record.qname
                break

def _subset_sequences( sequences, selected ):
    """
    Subset only the sequences that match 
    """
    for record in sequences:
        name = record.name.split()[0]
        if name in selected:
            yield record

def _write_sequences( sequences, output_file ):
    """
    Does what is says on the tin
    """
    with FastaWriter( output_file ) as writer:
        for record in sequences:
            writer.writeRecord( record )

if __name__ == '__main__':
    import sys

    input_file = sys.argv[1]
    reference_file = sys.argv[2]
    output_file = sys.argv[3]
    
    top_amp_analysis( input_file, reference_file, output_file )
