#! /usr/bin/env python

import sys

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbhla.io.BlasrIO import BlasrReader
from pbhla.fasta.utils import reverse_complement, write_fasta

def update_orientation( fasta_file, blasr_file, output_file ):
    fasta_records = list(FastaReader( fasta_file ))
    reversed_seqs = identify_reversed_sequences( blasr_file )
    reversed_records = reverse_records( fasta_records, reversed_seqs )
    write_fasta( reversed_records, output_file )
        
def identify_reversed_sequences( blasr_file ):
    reversed_seqs = []
    for record in BlasrReader( blasr_file ):
        if record.qstrand != record.tstrand:
            reversed_seqs.append( record.qname )
    return set(reversed_seqs)

def reverse_records( fasta_records, reversed_seqs ):
    reversed_records = []
    for record in fasta_records:
        name = record.name.split()[0]
        if name in reversed_seqs:
            record = reverse_complement( record )
        reversed_records.append( record )
    return reversed_records

if __name__ == '__main__':
    fasta_file = sys.argv[1]
    blasr_file = sys.argv[2]
    output_file = sys.argv[3]
    update_orientation( fasta_file, blasr_file, output_file )
