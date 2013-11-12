#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

from pbcore.io.FastqIO import FastqReader, FastqWriter, FastqRecord

def rename_resequencing( input_fastq, output_fastq ):
    """
    Rename resequenced Fastq to have an AA-style NumReads tag
    """
    with FastqWriter( output_fastq ) as writer:
        for record in FastqReader( input_fastq ):
            new_name = get_new_name( record.name )
            new_record = FastqRecord(new_name,
                                     record.sequence,
                                     record.quality)
            writer.writeRecord( new_record )

def get_new_name( name ):
    """
    Move the num_frags tag into the name as NumReads
    """
    parts = name.split()
    num_frags_parts = [p for p in parts if p.startswith('num_frags')]
    assert len(num_frags_parts) == 1
    num_frags = int(num_frags_parts[0].split('=')[1])
    return parts[0] + '_NumReads' + str(5*num_frags)