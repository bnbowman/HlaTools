#! /usr/bin/env python

import os, logging

from pbcore.io.FastaIO import FastaReader
from pbhla.io.BlasrIO import BlasrReader
from pbhla.fasta.utils import extract_names
from pbhla.annotation.HlaType import HlaType

log = logging.getLogger()

def summarize_alleles( read_align, ref_align, white_list, output=None ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    if output is None:
        basename = '.'.join( read_align.split('.')[:-1] )
        output = '%s.summary' % basename
    counts = _count_read_alignments( read_align )
    loci = _parse_reference_loci( ref_align )
    loci = _add_artificial_read_loci( loci, counts )
    # Combine the counts in various ways
    locus_counts = _count_reads_by_locus( counts, loci )
    white_list_ids = set(extract_names( white_list ))
    chimera_counts = _count_chimeras_by_locus( counts, loci, white_list_ids )
    nonchimera_counts = _calculate_nonchimeras( locus_counts, chimera_counts )
    # Summarize the counts for output
    read_summaries = _summarize_white_listed( counts, loci, white_list_ids, 
                                                            nonchimera_counts, 
                                                            chimera_counts )
    locus_summaries = _summarize_loci( locus_counts, nonchimera_counts, 
                                                     chimera_counts )
    _write_count_summary( loci, read_summaries, locus_summaries, output )

def _count_read_alignments( read_align ):
    """
    Count the number of reads associated with
    """
    counts = {}
    for record in BlasrReader( read_align ):
        try:
            counts[record.tname] += 1
        except:
            counts[record.tname] = 1
    return counts

def _parse_reference_loci( ref_align ):
    """
    Identify the loci associated with each reference sequence
    """
    loci = {}
    for record in BlasrReader( ref_align ):
        typing = HlaType.from_string( record.tname )
        loci[record.qname] = typing.gene
    return loci

def _add_artificial_read_loci( loci, counts ):
    """
    Add the locus for each artificial read to the locus dictionary
    """
    for read in counts:
        if read.startswith('Locus'):
            prefix = read.split('_')[0]
            locus = prefix[5:]
            loci[read] = locus
        elif read.startswith('HLA'):
            prefix = read.split('*')[0]
            locus = prefix.split('_')[-1]
            loci[read] = locus
    return loci

def _count_reads_by_locus( counts, loci ):
    """
    Sum up the counts of reads associated with each loci
    """
    locus_counts = {}
    for read, count in counts.iteritems():
        locus = loci[read]
        try:
            locus_counts[locus] += count
        except:
            locus_counts[locus] = count
    return locus_counts

def _count_chimeras_by_locus( counts, loci, white_list ):
    """
    Sum up the counts of chimeras associated with each loci
    """
    chimera_counts = {}
    for read, count in counts.iteritems():
        if read in white_list:
            continue
        locus = loci[read]
        try:
            chimera_counts[locus] += count
        except:
            chimera_counts[locus] = count
    return chimera_counts

def _calculate_nonchimeras( locus_counts, chimera_counts ):
    """
    Calculate the number of non-chimeras from the Total and Chimera counts
    """
    nonchimera_counts = {}
    for locus, count in locus_counts.iteritems():
        if locus in chimera_counts:
            nonchimera_counts[locus] = count - chimera_counts[locus]
        else:
            nonchimera_counts[locus] = count
    return nonchimera_counts

def _summarize_white_listed( counts, loci, white_list_ids, nonchimeras, chimeras ):
    """
    Summarize the counts and percentages for each white-listed read
    """
    summaries = {}
    total = float(sum([v for k,v in counts.iteritems()]))
    nonchimera_total = float(sum([v for k,v in nonchimeras.iteritems()]))
    for read in white_list_ids:
        locus = loci[read]
        summaries[read] = { 'Count':          counts[read],
                            'Total':          counts[read] + chimeras.get(locus, 0),
                            'AlleleFrac':     round(100*counts[read]/float(nonchimeras[locus]), 2),
                            'NonChimeraFrac': round(100*counts[read]/nonchimera_total, 2),
                            'AbsoluteFrac':   round(100*counts[read]/total, 2) }
    return summaries

def _summarize_loci( locus_counts, nonchimeras, chimeras ):
    """
    Summarize the counts and percentages for each locus
    """
    summaries = {}
    total = float(sum([v for k,v in locus_counts.iteritems()]))
    nonchimera_total = float(sum([v for k,v in nonchimeras.iteritems()]))
    for locus, count in nonchimeras.iteritems():
        allele_total = float(locus_counts[locus])
        summaries[locus] = { 'Count':          count,
                             'Total':          int(allele_total),
                             'AlleleFrac':     round(100*count/allele_total, 2),
                             'NonChimeraFrac': round(100*count/nonchimera_total, 2),
                             'AbsoluteFrac':   round(100*count/total, 2) }
    summaries['All'] = { 'Count':          int(nonchimera_total),
                         'Total':          int(total),
                         'AlleleFrac':     'N/A',
                         'NonChimeraFrac': round(100*nonchimera_total/total, 2),
                         'AbsoluteFrac':   100.0 }
    return summaries

def _write_count_summary( loci, read_summaries, locus_summaries, output ):
    with open( output, 'w' ) as handle:
        handle.write('Sequence\tCount\tTotal\tAlleleFrac\tNonChimeraFrac\tAbsoluteFrac\n')
        for locus in ['A', 'B', 'C', 'H', 'All']:
            for read, summary in read_summaries.iteritems():
                if loci[read] == locus:
                    handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (read,
                                                               summary['Count'],
                                                               summary['Total'],
                                                               summary['AlleleFrac'],
                                                               summary['NonChimeraFrac'],
                                                               summary['AbsoluteFrac']))
            if locus in locus_summaries:
                summary = locus_summaries[locus]
                handle.write('%s\t\t\t\t\t%s\t%s\t%s\t%s\t%s\n' % (locus,
                                                                   summary['Count'],
                                                                   summary['Total'],
                                                                   summary['AlleleFrac'],
                                                                   summary['NonChimeraFrac'],
                                                                   summary['AbsoluteFrac']))


if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    read_align = sys.argv[1]
    ref_align = sys.argv[2]
    white_list = sys.argv[3]
    
    summarize_alleles( read_align, ref_align, white_list )
