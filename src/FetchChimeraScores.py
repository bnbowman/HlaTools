#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
from operator import itemgetter

from pbcore.io import FastaReader, FastaWriter
from pbhla.utils import get_barcode, create_directory, check_output_file
from pbhla.external.utils import full_align_best_reference
from pbhla.io.BlasrIO import BlasrReader, pctsimilarity
from pbhla.fasta.utils import fasta_size

def fetch_chimera_scores( folder, reference ):
    # Check and read the reference data
    reference_files = find_reference_files( reference )

    # Check and read the query data
    barcode_files = split_results( folder )

    # Check and read the log data
    log_data = read_log_data( folder )

    good = []
    bad = []
    for barcode, reference in reference_files.iteritems():
        sample_file = barcode_files[barcode]
        alignments = full_align_best_reference( sample_file, reference )
        by_reference = hits_by_reference( alignments )
        best, rest = separate_best_hits( by_reference )
        good += best
        bad += rest

    print "Found %s 'good' and %s 'bad' consensus sequences" % (len(good), len(bad))
    good_scores = sorted({r.qname: log_data[r.qname] for r in good}.iteritems(), key=itemgetter(1), reverse=True)
    bad_scores = sorted({r.qname: log_data[r.qname] for r in bad}.iteritems(), key=itemgetter(1))
    print [k for k in good_scores[:5]]
    print [k for k in bad_scores]

def hits_by_reference( alignment ):
    by_reference = {}
    for hit in BlasrReader(alignment):
        try:
            by_reference[hit.tname].append( hit )
        except:
            by_reference[hit.tname] = [ hit ]
    return by_reference

def separate_best_hits( by_reference ):
    best, rest = [], []
    for reference, hits in by_reference.iteritems():
        if len(hits) == 1:
            best.append(hits[0])
        else:
            by_pctid = sorted(hits, key=pctsimilarity, reverse=True)
            best.append(by_pctid[0])
            rest += by_pctid[1:]
    return (best, rest)


def split_results( amp_analysis ):
    """Split the output of an Amplicon Analysis job by Barcode"""
    assert os.path.isdir( amp_analysis )
    sequence_path = os.path.join( amp_analysis, 'amplicon_analysis.fasta')
    check_output_file( sequence_path )
    print "Analyzing %s output sequences" % fasta_size( sequence_path )
    barcode_path = os.path.join( amp_analysis, 'by_barcode' )
    create_directory( barcode_path )

    records = list(FastaReader( sequence_path ))
    barcodes = {get_barcode(r):[] for r in records}
    [barcodes[get_barcode(r)].append(r) for r in records]
    barcode_files = {}
    for barcode, records in barcodes.iteritems():
        barcode_file = barcode + '.fasta'
        sample_path = os.path.join( barcode_path, barcode_file )
        with FastaWriter( sample_path ) as handle:
            for record in records:
                handle.writeRecord( record )
        barcode_files[barcode] = sample_path
    return barcode_files

def find_reference_files( reference_path ):
    assert os.path.isdir( reference_path )
    reference_files = {}
    for filename in os.listdir( reference_path ):
        if '--' in filename:
            barcode = filename.split('.')[0]
            reference_files[barcode] = os.path.join( reference_path, filename )
    return reference_files

def read_log_data( folder ):
    chimera_scores = {}
    log_path = os.path.join( folder, 'amplicon_analysis.log')
    check_output_file( log_path )

    with open( log_path ) as handle:
        for line in handle:
            line_parts = line.strip().split()
            consensus = line_parts[2]
            if 'abundant,' in line_parts:
                chimera_scores[consensus] = 0.0
            elif 'chimera' in line_parts:
                chimera_scores[consensus] = float(line_parts[8])
    return chimera_scores



if __name__ == '__main__':
    import sys

    folder = sys.argv[1]
    reference = sys.argv[2]

    fetch_chimera_scores( folder, reference )