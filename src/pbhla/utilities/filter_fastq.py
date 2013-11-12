#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

import logging
from pbcore.io.FastqIO import FastqWriter, FastqReader

log = logging.getLogger(__name__)

def filter_fastq( input_fastq, output_fastq, min_length=None, min_num_reads=None ):
    """
    Filter a Fastq file based on various criteria
    """
    kept = 0
    total = 0
    with FastqWriter( output_fastq ) as writer:
        for record in FastqReader( input_fastq ):
            total += 1
            if min_length and len(record.sequence) < min_length:
                continue
            if min_num_reads and get_num_reads( record ) < min_num_reads:
                continue
            kept += 1
            writer.writeRecord( record )
    log.info("Kept %s of %s consensus sequences" % (kept, total))

def get_num_reads( record ):
    return int( record.name.split('NumReads')[-1] )