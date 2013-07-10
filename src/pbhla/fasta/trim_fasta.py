#! /usr/bin/env python

import logging

from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord
from pbhla.io.BlasrIO import BlasrReader
from pbhla.fasta.utils import write_fasta

log = logging.getLogger(__name__)

WINDOW = 100
LOCI = ['A', 'B', 'C', 'DPA', 'DQA', 'DQB']

def trim_fasta( fasta_file, blasr_file, output_file, locus_dict, window=WINDOW, loci=LOCI ):
    log.info('Trimming sequences in "%s"' % fasta_file)
    log.debug("\tWindow Size:\t%s" % window)

    records = list( FastaReader( fasta_file ) )
    trims = parse_trims( blasr_file, window )
    trims = filter_trims_on_loci( trims, locus_dict, loci )
    trimmed_records = apply_trims( records, trims )
    write_fasta( trimmed_records, output_file )

    log.info('Finished trimming the supplied sequencs\n')
    return 

def parse_trims( blasr_file, window ):
    trims = {}
    for record in BlasrReader( blasr_file ):
        start = max(int(record.qstart)-window, 0)
        end = min(int(record.qend)+window, int(record.qlength))
        trims[record.qname] = (start, end)
    return trims

def filter_trims_on_loci( trims, locus_dict, loci ):
    filtered_trims = {}
    for key in trims:
        if locus_dict[key] in loci:
            filtered_trims[key] = trims[key]
    return filtered_trims

def apply_trims( records, trims ):
    trimmed = []
    for record in records:
        name = record.name.split()[0]
        if name in trims:
            start, end = trims[name]
            trimmed_record = FastaRecord( name, record.sequence[start:end] )
            trimmed.append( trimmed_record )
        else:
            trimmed.append( record )
    return trimmed
