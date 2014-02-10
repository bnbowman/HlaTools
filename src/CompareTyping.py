#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
import re
from collections import Counter

def compare_typing( typing_files ):
    files = sort_files( typing_files )
    barcodes = parse_barcodes( files )
    print "Indentified %s Barcodes to analyze" % len(barcodes)

    for bc in barcodes:
        print
        print 'Analyzing Barcode "%s"...' % bc
        records = extract_records( files, bc )

        if test_counts( records ) > 1:
            print "Uneven results detected for %s, skipping..." % bc
            continue
        else:
            print "Even results detected for all replicates of %s" % bc

        if test_genomic( records ) > 1:
            print "Differing genomic references detected for %s" % bc
            continue
        else:
            print "Matching genomic references found for all replicates of %s" % bc

        if test_cDNA( records ) > 1:
            print "Differing cDNA references detected for %s" % bc
            continue
        else:
            print "Matching cDNA references found for all replicates of %s" % bc

        if test_type( records ) > 1:
            print "Differing final-typings detected for %s" % bc
            continue
        else:
            print "Matching final-typings found for all replicates of %s" % bc

        if test_final( records ) > 1:
            print "Differing records detected for %s" % bc
            continue
        else:
            print "Matching records references found for all replicates of %s" % bc

def sort_files( typing_files ):
    files = {}
    for filename in typing_files:
        filepath = os.path.abspath( filename )
        parts = filepath.split('/')
        if len(parts) > 1:
            name = parts[-2]
            assert name not in files
            files[name] = filepath
    return files

def parse_barcodes( files ):
    barcodes = set()
    for name, filepath in files.iteritems():
        with open( filepath ) as handle:
            for line in handle:
                if line.startswith('Sequence'):
                    continue
                if not line.strip():
                    continue
                barcode = line.split('_')[0][7:]
                barcodes.add( barcode )
    return sorted(barcodes, key=lambda x: int(x.split('R')[1]))

def extract_records( files, bc ):
    records = {}
    for name, filepath in files.iteritems():
        records[name] = []
        with open( filepath ) as handle:
            for line in handle:
                if re.search(bc, line):
                    records[name].append( line.strip() )
    return records

def test_counts( records ):
    counts = {name: len(records[name]) for name in records}
    counter = Counter(counts.itervalues())
    return len(counter)

def test_genomic( records ):
    genomic_refs = {name: tuple([r.split()[2] for r in records[name]]) for name in records}
    counter = Counter(genomic_refs.itervalues())
    return len(counter)

def test_cDNA( records ):
    genomic_refs = {name: tuple(sorted([r.split()[7] for r in records[name]])) for name in records}
    counter = Counter(genomic_refs.itervalues())
    return len(counter)

def test_type( records ):
    genomic_refs = {name: tuple(sorted([r.split()[8] for r in records[name]])) for name in records}
    counter = Counter(genomic_refs.itervalues())
    return len(counter)

def test_final( records ):
    genomic_refs = {name: tuple(sorted([' '.join(r.split()[2:]) for r in records[name]])) for name in records}
    counter = Counter(genomic_refs.itervalues())
    #for c in counter:
    #    print c
    return len(counter)

if __name__ == '__main__':
    import sys

    typing_files = sys.argv[1:]

    compare_typing( typing_files )