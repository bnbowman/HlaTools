#! /usr/bin/env python

import logging, os

from pbhla.log import initialize_logger
from pbhla.references.data import get_exon_reference, get_genomic_reference, get_cDNA_reference
from pbhla.external.utils import full_align_best_reference
from pbhla.external.align_by_identity import align_by_identity
from pbhla.sequences.orientation import orient_sequences
from pbhla.sequences.input import get_input_file
from pbhla.sequences.rename import rename_sequences
from pbhla.alleles.extract import extract_alleles
from pbhla.cdna.extract_cDNA import extract_cDNA
from pbhla.typing.summarize import summarize_typing

GROUPINGS = ['locus', 'allele', 'both', 'all']
GROUPING = 'both'

log = logging.getLogger()

def type_sequences( input, grouping=GROUPING,
                           exon_fofn=None,
                           genomic_reference=None,
                           cDNA_reference=None,
                           loci=None):
    """
    Pick the top Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    log_file = get_log_file( input )
    initialize_logger( log, log_file=log_file )

    # First, get any references not specified by the user
    grouping = grouping or GROUPING
    exon_fofn = exon_fofn or get_exon_reference()
    genomic_reference = genomic_reference or get_genomic_reference()
    cDNA_reference = cDNA_reference or get_cDNA_reference()

    # Second, get the input file if a directory was specified
    sequence_file = get_input_file( input )

    # Finally, run the Typing procedure
    renamed_file = rename_sequences( sequence_file )
    raw_alignment = full_align_best_reference( renamed_file, genomic_reference )
    reoriented = orient_sequences( renamed_file, alignment_file=raw_alignment )
    selected = extract_alleles( reoriented, alignment_file=raw_alignment,
                                            method=grouping,
                                            loci=loci)
    gDNA_alignment = full_align_best_reference( selected, genomic_reference )
    cDNA_file = extract_cDNA( selected, exon_fofn, alignment_file=gDNA_alignment )
    cDNA_alignment = align_by_identity( cDNA_file, cDNA_reference )
    typing = summarize_typing( gDNA_alignment, cDNA_alignment )
    return typing

def get_output_dir( input ):
    if os.path.isdir( input ):
        return input
    else:
        return os.path.split( input )[0]

def get_log_file( input ):
    dir = get_output_dir( input )
    return os.path.join(dir, 'amp_analysis_typing.log')