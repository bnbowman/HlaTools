#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

import logging

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbhla.io.BlasrIO import BlasrReader

log = logging.getLogger(__name__)

def separate_amplicons( subread_file, alignment_file ):
    """
    Stuff
    """
    locations = _parse_alignment( alignment_file )
    medians = _calculate_medians( locations )
    centroids = _identify_centroids( locations, medians )
    assignments = _assign_reads( medians, centroids )
    _write_assigned_reads( subread_file, assignments )

def _parse_alignment( alignment ):
    """
    Parse the location of each hit in the alignment file
    """
    log.info("Parsing subread locations from alignment data")
    locations = {}
    for entry in BlasrReader( alignment ):
        if entry.tstrand == '1':
            start = int(entry.tlength) - int(entry.tend)
            end = int(entry.tlength) - int(entry.tstart)
        else:
            start = int(entry.tstart)
            end = int(entry.tend)
        locations[entry.qname] = (start, end)
    return locations

def _calculate_medians( locations ):
    """
    Convert a Dict of Start/End pairs to a Dict of Medians
    """
    return {k: (l[0]+l[1])/2 for k, l in locations.iteritems()}

def _identify_centroids( locations, medians ):
    """
    Identify the center of each of the two amplicons
    """
    log.info("Identifying the centroid of each amplicon")
    min_pos = min([s for s, e in locations.itervalues()])
    max_pos = max([e for s, e in locations.itervalues()])
    mid_pos = (min_pos + max_pos) / 2
    five_prime, three_prime = _split_medians( medians, mid_pos )
    #five_prime_center = _calculate_centroid( five_prime )
    #three_prime_center = _calculate_centroid( three_prime )
    five_prime_center = (min_pos + mid_pos) / 2
    three_prime_center = (max_pos + mid_pos) / 2
    return (five_prime_center, three_prime_center)

def _split_medians( medians, cutoff ):
    """
    Split the 5' and 3' reads
    """
    five_prime = []
    three_prime = []
    for key, median in medians.iteritems():
        if median < cutoff:
            five_prime.append( median )
        elif median > cutoff:
            three_prime.append( median )
    return five_prime, three_prime

def _calculate_centroid( positions ):
    """
    Calculate the 'Peak' or modal position from a set of reads
    """
    binned_pos = [p/BIN_SIZE*BIN_SIZE for p in positions]
    min_bin = -1 * BIN_SIZE
    max_bin = max(positions) + 2 * BIN_SIZE
    bins = range(min_bin, max_bin, BIN_SIZE)
    counts = [binned_pos.count(b) for b in bins]
    scores = [counts[c]-max(counts[c-1], counts[c+1])
                       for c in range(1,len(counts)-1)]
    sort_scores = sorted( scores, reverse=True )
    peak_bin_pos = [scores.index(s)+1 for s in sort_scores[:1]]
    peak_bins = [bins[p] for p in peak_bin_pos]
    peak_ranges = [(p-BIN_SIZE, p+2*BIN_SIZE) for p in peak_bins]
    peak_pos = [subset_positions(positions, start, end) for start, end in peak_ranges]
    peak_loc = [median(p) for p in peak_pos]
    return peak_loc[0]

def _assign_reads( medians, centroids ):
    """
    Assign reads to the centroid they are closer to
    """
    log.info("Assigning subreads reads to the closet amplicon cluster")
    assignments = {'5p':set(), '3p':set()}
    five_prime, three_prime = centroids
    for read, median in medians.iteritems():
        five_prime_diff = abs(median - five_prime)
        three_prime_diff = abs(median - three_prime)
        if five_prime_diff < three_prime_diff:
            assignments['5p'].add( read )
        else:
            assignments['3p'].add( read )
    return assignments


def _write_assigned_reads( input_fasta, assignments ):
    """
    Write out subreads to the appropriate file
    """
    log.info("Separating subreads based on their amplicon assignments")
    output_files = []
    writers = {}
    root_name = '.'.join( input_fasta.split('.')[:-1] )
    # Open up output writers for each group
    for group in assignments:
        output_file = "%s_%s.fasta" % (root_name, group)
        output_files.append( output_file )
        writers[group] = FastaWriter( output_file )
    # Write each record to it's appropriate group(s)
    for record in FastaReader( input_fasta ):
        name = record.name.split()[0]
        for group in assignments:
            if name in assignments[group]:
                writers[group].writeRecord( record )
                break
    # Close all of the output writers
    for group in writers:
        writers[group].close()
    return output_files

if __name__ == '__main__':
    import sys

    subread_file = sys.argv[1]
    alignment_file = sys.argv[2]

    logging.basicConfig( stream=sys.stdout )
    separate_amplicons( subread_file, alignment_file )