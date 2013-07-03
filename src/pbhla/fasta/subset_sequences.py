#! /usr/bin/env python

import sys

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbhla.fasta.utils import write_fasta

def subset_sequences( fasta_file, summary_file, output_file ):
    seq_ids = identify_sequences( summary_file )
    sequences = subset_sequence_records( fasta_file, seq_ids )
    write_fasta( sequences, output_file )

def identify_sequences( summary_file ):
    sequence_ids = []
    for line in open( summary_file ):
        if line.startswith('Locus'):
            continue
        sequence_ids.append( line.split()[1] )
    return set(sequence_ids)

def subset_sequence_records( fasta_file, seq_ids ):
    sequences = []
    for record in FastaReader( fasta_file ):
        name = record.name.split()[0]
        name = name.split('|')[0]
        if name.endswith('_cns'):
            name = name[:-4]
        if name in seq_ids:
            sequences.append( record )
    return sequences
        
if __name__ == '__main__':
    fasta_file = sys.argv[1]
    summary_file = sys.argv[2]
    output_file = sys.argv[3]
    subset_sequences( fasta_file, summary_file, output_file )
