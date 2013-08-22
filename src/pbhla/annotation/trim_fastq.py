#! /usr/bin/env python

import sys
from pbcore.io.FastqIO import FastqReader, FastqWriter, FastqRecord

WINDOW = 10

def trim_fastq( fastq_file, output_file, window=WINDOW ):
    with FastqWriter( output_file ) as writer:
        for record in FastqReader( fastq_file ):
            start = _find_start( record, window )
            end = _find_end( record, window )
            trimmed_record = _trim_fastq( record, start, end )
            writer.writeRecord( trimmed_record )

def _find_start( record, window ):
    for i in range( window ):
        if record.quality[window-i] == 0:
            return window-i+1

def _find_end( record, window ):
    length = len(record.sequence)
    for i in range( window ):
        if record.quality[length-window+i] == 0:
            return length-window+i

def _trim_fastq( record, start, end ):
    trimmed_sequence = record.sequence[start:end]
    trimmed_quality = record.quality[start:end]
    return FastqRecord( record.name,
                        trimmed_sequence,
                        trimmed_quality )

if __name__ == '__main__':
    fastq_file = sys.argv[1]
    output_file = sys.argv[2]

    trim_fastq( fastq_file, output_file )
