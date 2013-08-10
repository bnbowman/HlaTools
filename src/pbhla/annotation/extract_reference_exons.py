#! /usr/bin/env python

import sys, csv

from collections import namedtuple
from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter

AfaInfo = namedtuple('AfaInfo', 'pos0 pos1 region codon')

def extract_exons( afa_file, info_file):
    records = list(_parse_aligned_fasta( afa_file ))
    regions = list(_parse_info_file( info_file ))
    for exon in _select_exons( regions ):
        output_file = _get_output_file( info_file, exon )
        with FastaWriter( output_file ) as output:
            for record in _extract_fasta_region( records, exon ):
                if len(set(record.sequence)) == 1:
                    continue
                output.writeRecord( record )
    
def _parse_aligned_fasta( afa_file ):
    for i, record in enumerate(FastaReader( afa_file )):
        if i == 0:
            length = len(record.sequence)
        if len(record.sequence) == length:
            yield record
        else:
            raise ValueError( 'Uneven sequence length!')

def _parse_info_file( info_file ):
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

def _select_exons( regions ):
    for region in regions:
        if region[0].lower().startswith('exon'):
            yield region

def _get_output_file( info_file, region ):
    basename = info_file[:-5]
    region_name = region[0]
    return '%s_%s.fasta' % (basename, region_name)

def _extract_fasta_region( records, region ):
    name, start, end = region
    print name, start, end
    for record in records:
        yield FastaRecord( record.name,
                           record.sequence[start:end] )

if __name__ == '__main__':
    afa_file = sys.argv[1]
    info_file = sys.argv[2]
    extract_exons( afa_file, info_file )
