#! /usr/bin/env python

import os
import logging

from pbhla.external.utils import full_align_best_reference
from pbhla.external.align_by_identity import align_by_identity
from pbhla.sequences.orientation import orient_sequences
from pbhla.cdna.extract_cDNA import extract_cDNA
from pbhla.typing.summarize import summarize_typing

log = logging.getLogger()

def type_sequences( input_file, exon_fofn, genomic_reference, cDNA_reference ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    reoriented = orient_sequences( input_file, genomic_reference )
    gDNA_alignment = full_align_best_reference( reoriented, genomic_reference )
    cDNA_file = extract_cDNA( reoriented, exon_fofn, alignment_file=gDNA_alignment )
    cDNA_alignment = align_by_identity( cDNA_file, cDNA_reference )
    summarize_typing( gDNA_alignment, cDNA_alignment )

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_file = sys.argv[1]
    exon_fofn = sys.argv[2]
    genomic_reference = sys.argv[3]
    cDNA_reference = sys.argv[4]
    
    type_sequences( input_file, exon_fofn, genomic_reference, cDNA_reference )
