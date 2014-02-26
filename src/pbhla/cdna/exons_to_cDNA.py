#! /usr/bin/env python

import sys, os, logging

from pbcore.io.FastaIO import FastaReader, FastaRecord
from pbcore.io.FastqIO import FastqReader, FastqRecord
from pbhla.filenames import get_file_type
from pbhla.sequences.utils import write_sequences

log = logging.getLogger()

def exons_to_cDNA( exon_file ):
    """
    Combine a multi-Fasta of Exon sequences into a mock cDNA
    """
    output_type = get_file_type( exon_file )
    output_file = _get_output_file( exon_file, output_type )
    records = _parse_exon_records( exon_file, output_type )
    log.info("Combinging %s exons sequences to cDNA" % len(records))
    if len( records ):
        sorted_records = _sort_records( records )
        cDNA_record = _combine_records( sorted_records )
        log.info("Writing cDNA sequence out to %s" % output_file)
        write_sequences( cDNA_record, output_file )

def _parse_exon_records( exon_file, output_type ):
    if output_type == 'fasta':
        return list( FastaReader( exon_file ))
    elif output_type == 'fastq':
        return list( FastqReader( exon_file ))
    msg = 'Exon data must be in either Fasta or Fastq format'
    log.error( msg )
    raise TypeError( msg )

def _get_output_file( exon_file, output_type ):
    folder, filename = os.path.split( exon_file )
    return os.path.join( folder, 'cDNA.%s' % output_type)

def _sort_records( fasta_records ):
    """
    Sort a list of Fasta records by ascending Exon number
    """
    return sorted( fasta_records, key=get_exon_num )

def _combine_records( records ):
    """
    Combine an order series of Exon records in to a cDNA record
    """
    name = '_'.join( records[0].name.split('_')[:-1] )
    cDNA_sequence = ''
    cDNA_quality = ''
    for record in records:
        cDNA_sequence += record.sequence
        if hasattr( record, 'qualityString' ):
            cDNA_quality += record.qualityString
    if len(cDNA_sequence) == len(cDNA_quality):
        return FastqRecord( name, cDNA_sequence, qualityString=cDNA_quality )
    else:
        return FastaRecord( name, cDNA_sequence )

def get_exon_num( record ):
    return int(record.name[-1])

if __name__ == '__main__':
    exon_fasta = sys.argv[1]

    logging.basicConfig( stream=sys.stdout, level=logging.INFO )

    exons_to_cDNA( exon_fasta )
