#! /usr/bin/env python

import sys, os, logging

from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter
from pbhla.fasta.utils import write_fasta

log = logging.getLogger()

def exons_to_cDNA( exon_fasta ):
    """
    Combine a multi-Fasta of Exon sequences into a mock cDNA
    """
    output_file = _get_output_file( exon_fasta )
    fasta_records = list( FastaReader( exon_fasta ))
    if len(fasta_records):
        sorted_records = _sort_fasta_records( fasta_records )
        cDNA_record = _combine_records( sorted_records )
        write_fasta( [cDNA_record], output_file )

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
    exon_fasta = sys.argv[1]

    logging.basicConfig( stream=sys.stdout, level=logging.INFO )

    exons_to_cDNA( exon_fasta )
