#! /usr/bin/env python

import os, logging
from operator import itemgetter

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbcore.io.FastqIO import FastqReader, FastqWriter
from pbhla.fasta.utils import fasta_size, write_fasta
from pbhla.external.utils import get_alignment_file
from pbhla.io.BlasrIO import BlasrReader
from pbhla.utils import get_file_type, check_output_file

NPROC = 6
LOCI = ['A', 'B', 'C']
log = logging.getLogger()

def extract_alleles( input_file, output_file=None, reference_file=None, alignment_file=None, loci=LOCI ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    alignment_file = get_alignment_file( input_file, reference_file, alignment_file )
    # Set the output file if not specified
    output_file = output_file or _get_output_file( input_file )
    output_type = get_file_type( output_file )
    # Parse the alignment data and extract the target sequences
    alignments = list( BlasrReader( alignment_file ))
    groups = _group_by_locus( alignments, loci )
    ordered = _sort_groups( groups )
    selected = list( _select_sequences( ordered ))
    sequences = _parse_input_records( input_file )
    subset = list( _subset_sequences( sequences, selected ))
    _write_output( subset, output_file, output_type )
    return output_file

def _group_by_locus( alignments, loci ):
    """
    Group reads by the locus of their best alignment
    """
    groups = {}
    for record in alignments:
        reference = record.tname.split('*')[0]
        locus = reference.split('_')[-1]
        if locus not in loci:
            continue
        try:
            groups[locus].append( record )
        except:
            groups[locus] = [record]
    return groups

def _sort_groups( groups ):
    """
    Order each group of records individually
    """
    ordered = {}
    for locus, group in groups.iteritems():
        ordered[locus] = _sort_group( group )
    return ordered

def _sort_group( group ):
    """
    Order records in a group by their number of reads
    """
    # Internal function for simplicity
    def record_size( record ):
        return int( record.qname.split('NumReads')[-1] )

    # Count, sort and return
    counts = {record: record_size(record) for record in group}
    tuples = sorted( counts.iteritems(), key=itemgetter(1), reverse=True )
    return [t[0] for t in tuples]

def _select_sequences( groups ):
    """
    Select the top 1-2 sequences for each Locus
    """
    for group in groups.itervalues():
        first, the_rest = group[0], group[1:]
        # Yeild the first sequence from each group
        yield first.qname
        # If there are only two, yield the second as well
        if len(the_rest) == 1:
            yield the_rest[0].qname
        # If there are more than two, yield the next sequence with a different reference
        else:
            for record in the_rest:
                if record.tname != first.tname:
                    yield record.qname
                    break

def _parse_input_records( input_file ):
    """
    Parse the input sequence records with the appropriate pbcore Reader
    """
    input_type = get_file_type( input_file )
    if input_type == 'fasta':
        return list( FastaReader( input_file ))
    elif input_type == 'fastq':
        return list( FastqReader( input_file ))
    else:
        msg = 'Input file must be either Fasta or Fastq'
        log.error( msg )
        raise TypeError( msg )

def _subset_sequences( sequences, selected ):
    """
    Subset only the sequences that match 
    """
    for record in sequences:
        name = record.name.split()[0]
        if name in selected:
            yield record

def _write_output( records, output_file, output_type ):
    """
    Write the records out to file
    """
    if output_type == 'fasta':
        write_fasta( records, output_file )
    else:
        with FastqWriter( output_file ) as writer:
            for record in records:
                writer.writeRecord( record )
        check_output_file( output_file )

def _get_output_file( input_file ):
    basename = '.'.join( input_file.split('.')[:-1] )
    file_type = get_file_type( input_file )
    return '%s.selected.%s' % (basename, file_type)

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reference_file = sys.argv[3]
    
    extract_alleles( input_file, output_file, reference_file )
