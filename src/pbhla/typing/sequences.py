#! /usr/bin/env python

import logging, logging.config

from pbhla import __LOG__
from pbhla.external.utils import (align_best_reference, 
                                  full_align_best_reference)
from pbhla.external.align_by_identity import align_by_identity
from pbhla.sequences.orientation import orient_sequences
from pbhla.alleles.extract import extract_alleles
from pbhla.cdna.extract_cDNA import extract_cDNA
from pbhla.typing.summarize import summarize_typing

GROUPINGS = ['locus', 'allele', 'both', 'all']
GROUPING = 'both'

logging.config.fileConfig( __LOG__ )
log = logging.getLogger()

def type_sequences( sequence_file, exon_fofn, genomic_reference, cDNA_reference, grouping=GROUPING ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    # First we align the sequences to the reference and annotate typing
    raw_alignment = align_best_reference( sequence_file, genomic_reference )
    reoriented = orient_sequences( sequence_file, alignment_file=raw_alignment )
    selected = extract_alleles( reoriented, alignment_file=raw_alignment, method=grouping )
    gDNA_alignment = full_align_best_reference( selected, genomic_reference )
    cDNA_file = extract_cDNA( selected, exon_fofn, alignment_file=gDNA_alignment )
    cDNA_alignment = align_by_identity( cDNA_file, cDNA_reference )
    summarize_typing( gDNA_alignment, cDNA_alignment )

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_file = sys.argv[1]
    exon_fofn = sys.argv[2]
    genomic_reference = sys.argv[3]
    cDNA_reference = sys.argv[4]
    grouping = sys.argv[5] if len(sys.argv) > 5 else GROUPING

    assert grouping in GROUPINGS

    type_sequences( input_file, exon_fofn, genomic_reference, cDNA_reference, grouping )
