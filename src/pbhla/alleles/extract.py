#! /usr/bin/env python

import logging
from operator import itemgetter

from pbcore.io.FastaIO import FastaReader
from pbcore.io.FastqIO import FastqReader, FastqWriter
from pbhla.fasta.utils import write_fasta, get_num_reads
from pbhla.external.utils import get_alignment_file
from pbhla.io.BlasrIO import BlasrReader
from pbhla.utils import get_file_type, check_output_file

NPROC = 6
LOCI = ['A', 'B', 'C', 'DRB1', 'DQB1']
METHOD = 'locus'
MIN_FRAC = 0.15

log = logging.getLogger()

def extract_alleles( input_file, output_file=None, reference_file=None,
                                                   alignment_file=None,
                                                   method=METHOD,
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

    log.info('Selecting top sequences from %s according to the "%s" policy' % (input_file, method))
    ordered = _sort_groups( groups )
    selected = list( _select_sequences( ordered ))
    sequences = _parse_input_records( input_file )
    log.info('Selected %s sequences from %s total for further analysis' % (len(selected), len(sequences)))
    subset = list( _subset_sequences( sequences, selected ))
    log.info('Writing the selected sequences out to %s' % output_file)
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
        barcode = name.split('_')[0]
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

def _sort_groups( groups ):
    """Order each group of records individually"""
    ordered = {}
    for locus, group in groups.iteritems():
        ordered[locus] = _sort_group( group )
    return ordered

def _sort_group( group ):
    """Order records in a group by their number of reads"""
    # Internal function for simplicity
    def record_size( record ):
        return int( record.qname.split('NumReads')[-1] )

    # Count, sort and return
    counts = {record: record_size(record) for record in group}
    tuples = sorted( counts.iteritems(), key=itemgetter(1), reverse=True )
    return [t[0] for t in tuples]

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

def _parse_input_records( input_file ):
    """Parse the input sequence records with the appropriate pbcore Reader"""
    input_type = get_file_type( input_file )
    if input_type == 'fasta':
        return list( FastaReader( input_file ))
    elif input_type == 'fastq':
        return list( FastqReader( input_file ))
    else:
        msg = 'Input file must be either Fasta or Fastq'
        log.error( msg )
        raise TypeError( msg )

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
