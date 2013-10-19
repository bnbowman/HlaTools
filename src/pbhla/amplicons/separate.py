#! /usr/bin/env python

import re
import os
import logging

from pbcore.io.FastaIO import FastaReader, FastaWriter

from pbhla.external.commandline_tools import run_blasr
from pbhla.io.BlasrIO import BlasrReader
from pbhla.fasta.utils import fasta_size
from pbhla.utils import check_output_file, read_list_file, write_list_file
from pbhla.references.fofn import parse_reference_dict

BIN_SIZE = 100

log = logging.getLogger()

def separate_amplicons( input_data, reference_fofn, target_loci, output=None ):
    """
    Public interface for _separate_subreads
    """
    # Check and set the input and output, as needed
    log.info("Separating amplicons for these loci[%s]" % ','.join(target_loci))
    file_list = _parse_input( input_data )
    output = output or _get_output_file( input_data )
    references = parse_reference_dict( reference_fofn )

    # Iterate over the input subread files, splitting as needed
    new_files = []
    for filepath in file_list:
        locus = get_file_locus( filepath )
        if locus in target_loci:
            if is_amplicon_specific( filepath ):
                log.info("Subreads for Locus %s already split, skipping..." % locus)
                new_files.append( filepath )
            else:
                log.info("Subreads for Locus %s already split, skipping..." % locus)
                reference_fasta = references[locus]
                new_file_list = _separate_amplicons( file_list, reference_fasta, locus)
                continue

        # Otherwise, separate the sequences and write the results
        log.info("Separating subreads by amplicon for Locus %s" % locus)
        new_file_list = _separate_amplicons( file_list, reference_fasta, locus)
    write_list_file( new_file_list, output )

def _parse_input(input_data):
    """
    Parse the list of subread files from the input if needed
    """
    if isinstance(input_data, str):
        return read_list_file( input_data )
    elif isinstance(input_data, list):
        return input_data
    else:
        msg = 'Input must be FOFN or List'
        log.error( msg )
        raise TypeError( msg )

def _get_output_file(input_file):
    """
    If an output file wasn't specified, set default output file
    """
    if isinstance(input_file, str) and output is None:
        return input_file
    else:
        msg = 'Output file must be specified with file-list input!'
        log.error( msg )
        raise ValueError( msg )

def is_amplicon_specific( filepath ):
    dirname, basename = os.path.split( filepath )
    filename = basename.split('.')[0]
    if re.search('_3p_', filename) or re.search('_5p_', filename):
        return True
    return False

def is_allele_specific( filepath ):
    dirname, basename = os.path.split( filepath )
    filename = basename.split('.')[0]
    if re.search('_A1_', filename) or re.search('_A2_', filename):
        return True
    return False

def get_file_locus( filepath ):
    dirname, basename = os.path.split( filepath )
    filename = basename.split('.')[0]
    parts = filename.split('_')
    return parts[1]

def _separate_amplicons2( file_list, reference_fasta, locus):
    """
    Separate the 5' and 3' aligned subreads for a given locus
    """
    subread_file, other_files = _separate_file_list( file_list, locus )
    alignment = _align_subreads( subread_file, reference_fasta, locus )
    locations = _parse_alignment( alignment )
    os.remove( alignment )
    medians = _calculate_medians( locations )
    centroids = _identify_centroids( locations, medians )
    assignments = _assign_reads( medians, centroids )
    new_subread_files = _write_assigned_reads( subread_file, assignments )
    return new_subread_files + other_files

def _separate_file_list( file_list, target_locus ):
    """
    Separate the target fasta file from the FOFN of all subread files
    """
    log.info("Parsing locus-specific subread FOFN")
    target_fasta = None
    other_fasta = []
    print file_list, target_locus
    for filename in file_list:
        basename = filename.split('.')[0]
        locus = basename.split('_')[-1]
        if locus == target_locus and target_fasta is None:
            target_fasta = filename
        elif locus == target_locus:
            msg = 'Multiple files for target locus found!'
            log.error( msg )
            raise ValueError( msg )
        else:
            other_fasta.append( filename )
    if target_fasta is None:
        msg = 'No fasta file for target locus found!'
        log.error( msg )
        raise ValueError( msg )
    return ( target_fasta, other_fasta )

def _align_subreads( subread_fasta, reference_fasta, locus ):
    """
    Align all locus-specific subreads against the appropriate references
    """
    location = os.path.dirname( subread_fasta )
    alignment_file = os.path.join(location, 'temp.m1')
    subread_count = fasta_size( subread_fasta )
    reference_count = fasta_size( reference_fasta )
    blasr_args = {'nproc': 8,
                  'out': alignment_file,
                  'bestn': 1,
                  'nCandidates': reference_count,
                  'noSplitSubreads': True}
    log.info("Aligning %s reads against %s references for %s" % (subread_count, 
                                                                 reference_count,
                                                                 locus))
    run_blasr( subread_fasta, reference_fasta, blasr_args )
    check_output_file( alignment_file )
    return alignment_file

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

#
# Utilities
#
def subset_positions(positions, start, end):
    return [p for p in positions
                    if p >= start
                    if p <= end]

def median(positions):
    positions = sorted(positions)
    if len(positions) % 2:
        middle = len(positions)/2+1
    else:
        middle = len(positions)/2
    return positions[middle]
   


if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    subread_fofn = sys.argv[1]
    reference_fofn = sys.argv[2]
    locus = sys.argv[3]
    output_file = sys.argv[4]

    separate_amplicons( subread_fofn, reference_fofn, locus, output_file )
