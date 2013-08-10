#! /usr/bin/env python

import sys, os, logging

from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter
from pbhla.fasta.utils import write_fasta

log = logging.getLogger()

def pair_exon_files( fofn_file, overlap_exon ):
    """
    Pair the 5' and 3' amplicons of a gene base on 1 overlapping exon    
    """
    exon_files = list( _parse_fofn( fofn_file ))
    output_file = _get_output_file( exon_fasta )
    fasta_records = list( FastaReader( exon_fasta ))
    sorted_records = _sort_fasta_records( fasta_records )
    cDNA_record = _combine_records( sorted_records )
    write_fasta( [cDNA_record], output_file )

def _parse_fofn( fofn_file ):
    with open( fofn_file, 'r' ) as handle:
        for line in handle:
            yield line.strip()

def _get_output_file( exon_fasta ):
    basename = '_'.join( exon_fasta.split('_')[:-1] )
    return '%s_cDNA.fasta' % basename

def _sort_fasta_records( fasta_records ):
    """
    Sort a list of Fasta records by ascending Exon number
    """
    return sorted( fasta_records,
                   key=get_exon_num )

def _combine_records( records ):
    """
    Combine an order series of Exon records in to a cDNA record
    """
    basename = '_'.join( records[0].name.split('_')[:-1] )
    cDNA_name = '%s_cDNA' % basename
    cDNA_sequence = ''
    for record in records:
        cDNA_sequence += record.sequence
    return FastaRecord( cDNA_name, cDNA_sequence )

def get_exon_num( record ):
    return int(record.name[-1])

if __name__ == '__main__':
    fofn_file = sys.argv[1]
    overlap_exon = int( sys.argv[2] )

    logging.basicConfig( stream=sys.stdout, level=logging.INFO )

    pair_exon_files( fofn_file, overlap_exon )
