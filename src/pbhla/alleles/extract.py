#! /usr/bin/env python

import logging
from operator import itemgetter

from pbcore.io.FastqIO import FastqWriter
from pbhla.fasta.utils import write_fasta, get_num_reads
from pbhla.external.utils import get_alignment_file
from pbhla.filenames import get_file_type
from pbhla.io.BlasrIO import BlasrReader
from pbhla.utils import check_output_file
from pbhla.record import consensus_size, record_accuracy
from pbhla.sequences.utils import read_sequences

NPROC = 6
LOCI = ['A', 'B', 'C', 'DRB1', 'DQB1']
METHOD = 'locus'
SORT = 'accuracy'
MIN_FRAC = 0.15

log = logging.getLogger()

def extract_alleles( input_file, output_file=None, reference_file=None,
                                                   alignment_file=None,
                                                   method=METHOD,
                                                   sort=SORT,
                                                   loci=LOCI ):
    """Pick the top 2 Amplicon Analysis consensus seqs per group from a Fasta"""
    method = method or METHOD
    loci = loci or LOCI

    # Set the output file if not specified
    output_file = output_file or _get_output_file( input_file )
    output_type = get_file_type( output_file )

    # If align to reference for breaking ties
    alignment_file = get_alignment_file( input_file, reference_file, alignment_file )
    alignments = list( BlasrReader( alignment_file ))

    # Run the appropriate grouping
    if method == 'locus':
        groups = _group_by_locus( alignments, loci )
    elif method == 'barcode':
        groups = _group_by_barcode( alignments )
    elif method == 'both':
        groups = _group_by_both( alignments, loci )
    elif method == 'all':
        groups = {a.qname: [a] for a in alignments}
    else:
        msg = "Invalid Selection Metric: %s" % method
        log.error( msg )
        raise ValueError( msg )

    # Read the input sequences and use them to generate our sorting data
    sequences = read_sequences( input_file )
    if sort == 'num_reads':
        sorting_data = {s.name: consensus_size(s) for s in sequences}
    elif sort == 'accuracy':
        assert get_file_type(input_file) == 'fastq'
        sorting_data = {s.name: record_accuracy(s) for s in sequences}
    else:
        msg = "Invalid Sorting Metric: %s" % sort
        log.error( msg )
        raise ValueError( msg )

    log.info('Sorting sequences for selection according to "%s"' % sort)
    ordered = _sort_groups( groups, sorting_data )

    log.info('Selecting top sequences from %s according to the "%s" policy' % (input_file, method))
    selected = list( _select_sequences( ordered ))
    log.info('Selected %s sequences from %s total for further analysis' % (len(selected), len(sequences)))

    log.info('Writing the selected sequences out to %s' % output_file)
    subset = list( _subset_sequences( sequences, selected ))
    _write_output( subset, output_file, output_type )
    return output_file

def _group_by_locus( alignments, loci ):
    """Group reads by the locus of their best alignment"""
    groups = {}
    for record in alignments:
        reference = record.tname.split('*')[0]
        locus = reference.split('_')[-1]
        if locus not in loci:
            continue
        try:
            groups[locus].append( record )
        except KeyError:
            groups[locus] = [record]
    return groups

def _group_by_barcode( alignments ):
    """Group reads by their barcode"""
    groups = {}
    for alignment in alignments:
        name = alignment.qname
        if name.startswith('Barcode'):
            name = name[7:]
        if name.startswith('_'):
            name = name[1:]
        barcode = name.split('_Cluster')[0]
        try:
            groups[barcode].append( alignment )
        except KeyError:
            groups[barcode] = [ alignment ]
    return groups

def _group_by_both( alignments, loci ):
    groups = {}
    barcode_groups = _group_by_barcode( alignments )
    for barcode, bc_alignments in barcode_groups.iteritems():
        locus_groups = _group_by_locus( bc_alignments, loci )
        for locus, locus_alignments in locus_groups.iteritems():
            group = '%s_%s' % (barcode, locus)
            groups[group] = locus_alignments
    return groups

def _sort_groups( groups, sorting_data ):
    """Order each group of records individually"""
    ordered = {}
    for locus, group in groups.iteritems():
        print sorting_data
        sorted_records = sorted( group, key=lambda x: sorting_data[x.qname], reverse=True )
        ordered[locus] = sorted_records
    return ordered

def _select_sequences( groups ):
    """Select the top 1-2 sequences for each Locus"""
    for group in groups.itervalues():
        first, the_rest = group[0], group[1:]
        # Yield the first sequence from each group
        yield first.qname
        first_reads = get_num_reads( first.qname )
        # Yield the next sequence with a different reference
        for record in the_rest:
            num_reads = get_num_reads( record.qname )
            if record.tname == first.tname and record.nmis == first.nmis:
                first_reads += get_num_reads( record.qname )
            elif num_reads > (first_reads * MIN_FRAC):
                yield record.qname
                break

def _subset_sequences( sequences, selected ):
    """Subset only the sequences that match"""
    for record in sequences:
        name = record.name.split()[0]
        if name in selected:
            yield record

def _write_output( records, output_file, output_type ):
    """Write the records out to file"""
    if output_type == 'fasta':
        write_fasta( records, output_file )
    else:
        with FastqWriter( output_file ) as writer:
            for record in records:
                writer.writeRecord( record )
        check_output_file( output_file )

def _get_output_file( input_file ):
    basename = '.'.join( input_file.split('.')[:-1] )
    file_type = get_file_type( input_file )
    return '%s.selected.%s' % (basename, file_type)

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reference_file = sys.argv[3]
    
    extract_alleles( input_file, output_file, reference_file )
