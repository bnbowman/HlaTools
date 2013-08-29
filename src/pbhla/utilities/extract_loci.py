#! /usr/bin/env python

import os

from pbcore.io import FastaReader, FastaWriter

def extract_loci( project, loci, output ):
    contigs = list( _extract_contig_names( project, loci ))
    sequences = list( _extract_sequences( project, contigs ))
    _write_sequences( sequences, output )
    
def _extract_contig_names( project, loci ):
    summary_file = os.path.join( project, 'results', 'AmpliconAssembly', 'Locus_Calls.txt' )
    with open( summary_file, 'r' ) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            if line.startswith('Locus'):
                continue
            locus, contig, length, count, fraction, hit, pctid = line.strip().split()
            if locus in loci:
                yield contig

def _extract_sequences( project, contigs ):
    sequence_file = os.path.join( project, 'results', 'AmpliconAssembly', 'Final_Sequences.fasta' )
    for record in FastaReader( sequence_file ):
        name = record.name.split()[0]
        if name in contigs:
            yield record

def _write_sequences( sequences, output ):
    with FastaWriter( output ) as handle:
        for record in sequences:
            handle.writeRecord( record )

if __name__ == '__main__':
    import sys

    project_folder = sys.argv[1]
    loci = sys.argv[2].split(',')
    output = sys.argv[3]

    extract_loci( project_folder, loci, output )
