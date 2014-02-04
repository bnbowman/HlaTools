#! /usr/bin/env python

import os, logging

from pbhla.external.utils import (align_best_reference, 
                                  full_align_best_reference)
from pbhla.external.align_by_identity import align_by_identity
from pbhla.sequences.orientation import orient_sequences
from pbhla.io.AmpAnalysisUtils import (orient_amp_analysis,
                                       subset_amp_analysis)
from pbhla.alleles.extract import extract_alleles
from pbhla.cdna.extract_cDNA import extract_cDNA
from pbhla.typing.summarize import summarize_typing
from pbhla.typing.chimeras import create_chimeras

#from pbhla.io.extract_best_reads import extract_best_reads
#from pbhla.fasta.utils import combine_fasta
#from pbhla.annotation.summarize_alleles import summarize_alleles

log = logging.getLogger()

def type_sequences( input_folder, exon_fofn, genomic_reference, cDNA_reference ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    sequence_file = os.path.join( input_folder, 'amplicon_analysis.fastq' )
    csv_file = os.path.join( input_folder, 'amplicon_analysis.csv' )
    # First we align the sequences to the reference and annotate typing
    raw_alignment = align_best_reference( sequence_file, genomic_reference )
    reoriented = orient_sequences( sequence_file, alignment_file=raw_alignment )
    reoriented_csv = orient_amp_analysis( csv_file, raw_alignment )
    selected = extract_alleles( reoriented, alignment_file=raw_alignment )
    selected_csv = subset_amp_analysis( reoriented_csv, selected )
    gDNA_alignment = full_align_best_reference( selected, genomic_reference )
    cDNA_file = extract_cDNA( selected, exon_fofn, alignment_file=gDNA_alignment )
    cDNA_alignment = align_by_identity( cDNA_file, cDNA_reference )
    summarize_typing( gDNA_alignment, cDNA_alignment )
    # Next we generate some mock chimera sequences
    #chimera_file = create_chimeras( selected, alignment_file=gDNA_alignment )
    #basename = '.'.join( chimera_file.split('.')[:-2] )
    #combined_file = '%s.combined.fasta' % basename
    #combine_fasta( [input_file, chimera_file], combined_file )
    # Finally we use a competetive alignment of best-reads to summarize the allelic breakdown
    #dirname = os.path.dirname( input_fasta )
    #best_reads = os.path.join( dirname, 'reads_of_insert.fasta' )
    #extract_best_reads( input_fofn, best_reads )
    #best_alignment = align_best_reference( best_reads, combined_file )
    #summarize_alleles( best_alignment, raw_alignment, selected )

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_folder = sys.argv[1]
    exon_fofn = sys.argv[2]
    genomic_reference = sys.argv[3]
    cDNA_reference = sys.argv[4]
    
    type_sequences( input_folder, exon_fofn, genomic_reference, cDNA_reference )
