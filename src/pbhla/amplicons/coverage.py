#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

from pbhla.io.BlasrIO import BlasrReader

def calculate_coverage( alignment_file ):
    """
    Stuff
    """
    locations = _parse_alignment( alignment_file )
    coverage = {}
    for start, end, ref in locations.itervalues():
        if ref not in coverage:
            coverage[ref] = {}
        start = round_pos(start)
        for i in range(start, end, 10):
            try:
                coverage[ref][i] += 1
            except KeyError:
                coverage[ref][i] = 1
    for ref, data in coverage.iteritems():
        for key in sorted(data):
            print '%s\t%s\t%s' % (ref, str(key).ljust(5), data[key])


def _parse_alignment( alignment ):
    """
    Parse the location of each hit in the alignment file
    """
    locations = {}
    for entry in BlasrReader( alignment ):
        if entry.tstrand == '1':
            start = int(entry.tlength) - int(entry.tend)
            end = int(entry.tlength) - int(entry.tstart)
        else:
            start = int(entry.tstart)
            end = int(entry.tend)
        locations[entry.qname] = (start, end, entry.tname)
    return locations

def round_pos( pos ):
    return int(10 * round(pos/10.0, 0))

if __name__ == '__main__':
    import sys

    alignment_file = sys.argv[1]

    calculate_coverage( alignment_file )