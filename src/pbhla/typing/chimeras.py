#! /usr/bin/env python

import os, logging
from operator import itemgetter

from pbcore.io import FastaReader, FastaWriter, FastaRecord
from pbhla.fasta.utils import fasta_size, write_fasta
from pbhla.external.utils import align_best_reference
from pbhla.io.BlasrIO import BlasrReader

NPROC = 6

log = logging.getLogger()

def create_chimeras( input_file, output=None, reference_file=None, alignment_file=None ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    # Check the input files, and align the input file if needed
    if reference_file and alignment_file is None:
        alignment_file = align_best_reference( input_file, reference_file )
    elif reference_file is None and alignment_file is None:
        msg = "extract_alleles requires either an Alignment or a Reference!"
        log.error( msg )
        raise IOError( msg )
    # Set the output file if not specified
    if output is None:
        basename = '.'.join( input_file.split('.')[:-1] )
        output = '%s.chimeras.fasta' % basename
    # Parse the alignment data and extract the target sequences
    alignments = list( BlasrReader( alignment_file ))
    groups = _group_by_locus( alignments )
    groups = _filter_groups( groups )
    sequences = list( FastaReader( input_file ))
    chimeras = list( _create_chimeras( groups, sequences ))
    write_fasta( chimeras, output )
    return output

def _group_by_locus( alignments ):
    """
    Group reads by the locus of their best alignment
    """
    loci = {}
    for record in alignments:
        reference = record.tname.split('*')[0]
        locus = reference.split('_')[-1]
        try:
            loci[locus].append( record.qname )
        except:
            loci[locus] = [ record.qname ]
    return loci

def _filter_groups( groups ):
    """
    Filter homozygous loci, raise an error if 3+ alleles
    """
    filtered_groups = {}
    for locus, group in groups.iteritems():
        if len( group ) == 1:
            continue
        if len( group ) == 2:
            filtered_groups[locus] = group
        if len( group ) >= 3:
            msg = "Cannot create chimeras for more than 2 sequences"
            log.error( msg )
            raise ValueError( msg )
    return filtered_groups

def _create_chimeras( groups, sequences ):
    """
    Create chimeras for each homozygous locus group and return
    """
    chimeras = []
    for locus, group in groups.iteritems():
        subset = list( _subset_sequences( group, sequences ))
        chimeras += _create_chimera_pair( locus, subset )
    return chimeras

def _create_chimera_pair( locus, sequences ):
    """
    Create chimeras for a single locus
    """
    assert len( sequences ) == 2
    first_5p, first_3p = _split_sequence( sequences[0] )
    second_5p, second_3p = _split_sequence( sequences[1] )
    A_chimera = FastaRecord(name='Locus%s_ChimeraA' % locus,
                            sequence=first_5p+second_3p)
    B_chimera = FastaRecord(name='Locus%s_ChimeraB' % locus,
                            sequence=second_5p+first_3p)
    return [A_chimera, B_chimera]

def _split_sequence( record ):
    """
    Split a FastaRecord sequence in half and return
    """
    middle  = len( record.sequence ) / 2
    return (record.sequence[:middle], record.sequence[middle:])

def _subset_sequences( group, sequences ):
    """
    Subset a list of FastaRecords by name
    """
    for record in sequences:
        name = record.name.split()[0]
        if name in group:
            yield record

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_file = sys.argv[1]
    output = sys.argv[2]
    reference_file = sys.argv[3]
    
    create_chimeras( input_file, output, reference_file )
