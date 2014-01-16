#! /usr/bin/env python

import sys, csv

from collections import namedtuple
from pbhla.io.AlignIO import read_alignment, write_stockholm_alignment

AfaInfo = namedtuple('AfaInfo', 'pos0 pos1 region codon')

def extract_exons( afa_file, info_file, prefix=None ):
    """
    Extract all aligned regions corresponding to exons from an MSA
    """
    alignment = read_alignment( afa_file )
    prefix = prefix or '.'.join( afa_file.split('.')[:-1] )

    # Iterate over the Regions, extracting each Exon
    for region in _parse_regions( info_file ):
        print region
        if region[0].startswith('exon'):
            name, start, end = region
            output_file = '%s_%s.sto' % (prefix, name)
            exon_alignment = alignment[:,start:end]
            if len(exon_alignment):
                write_stockholm_alignment( output_file, exon_alignment )

def _parse_regions( info_file ):
    region = None
    start = None
    with open( info_file ) as handle:
        for record in map(AfaInfo._make, csv.reader( handle, delimiter=' ')):
            end = int( record.pos0 )
            if record.region != region and region:
                yield (region, start, end)
                start = int( record.pos0 )
            region = record.region
            start = int( record.pos0 ) if start is None else start
    yield (region, start, end)

if __name__ == '__main__':
    afa_file = sys.argv[1]
    info_file = sys.argv[2]
    extract_exons( afa_file, info_file )