#! /usr/bin/env python
import os, sys, shutil, logging

from pbcore.io.FastaIO import FastaReader 
from pbphase.clusense import Clusense

from pbhla.arguments import args, parse_args
from pbhla.fofn import create_baxh5_fofn
from pbhla.separate_sequences import (separate_sequences,
                                      separate_listed_sequences, 
                                      separate_aligned_sequences, 
                                      write_group,
                                      write_all_groups)
from pbhla.fasta.utils import (fasta_size, 
                               extract_sequence,
                               combine_fasta)
from pbhla.fasta.update_orientation import update_orientation
from pbhla.fasta.subset_sequences import subset_sequences
from pbhla.fasta.trim_fasta import trim_fasta
from pbhla.io.extract_subreads import extract_subreads
from pbhla.references import (create_fofn_reference,
                              create_m1_reference,
                              create_m5_reference,
                              create_sam_reference,
                              create_phased_reference,
                              filter_m5_file,
                              create_reference_fasta)
from pbhla.stats.SubreadStats import SubreadStats
from pbhla.phasing.SummarizeClusense import combine_clusense_output
from pbhla.resequencing.SummarizeResequencing import combine_resequencing_output 
from pbhla.annotation.summarize import ( summarize_contigs,
                                         meta_summarize_contigs,
                                         summarize_resequenced,
                                         summarize_typing )
from pbhla.annotation.typing import type_hla
from pbhla.external.HbarTools import HbarRunner
from pbhla.external.commandline_tools import run_blasr
from pbhla.external.CdHit import cd_hit_est
from pbhla.external.ClusterResequencer import ClusterResequencer
from pbhla.utils import *

log = logging.getLogger()

class NewSubreadPipeline( object ):

    def __init__( self ):
        # Parse the options 
        parse_args()
        # Initialize the rest of the project
        self._initialize_project()
        self._initialize_logging()

    def _initialize_project( self ):
        create_directory( args.output )
        # Create the various sub-directories
        for d in ['log_files', 'HBAR', 'references', 'subreads', 
                  'alignments', 'phasing', 'phasing_results', 
                  'resequencing', 'resequencing_results', 'results', 
                  'stats', 'annotation']:
            sub_dir = os.path.join( args.output, d )
            create_directory( sub_dir )
            setattr(self, d, sub_dir)

    def _initialize_logging( self ):
        log_format = "%(asctime)s [%(levelname)s] %(funcName)s %(message)s"
        log_file = os.path.join( self.log_files, "HLA_Pipeline.log" )
        logging.basicConfig( level=logging.INFO, 
                             format=log_format,
                             stream=sys.stdout )
        logging.basicConfig( level=logging.INFO, 
                             format=log_format, 
                             filename=log_file )
    
    def run( self ):
        self.extract_new_subreads()
        self.align_new_subreads_to_contigs()
        self.read_on_target_contigs()
        self.separate_off_target_subreads()
        self.create_ref_locus_dict()
        self.create_phased_locus_dict()
        self.align_new_subreads_to_ref()
        self.combine_subreads_by_allele()
        self.combine_subreads_by_locus()

    def extract_new_subreads( self ):
        log.info('Beginning the extraction of subread data')
        # Dump all valid reads from the above files
        self.subread_file = os.path.join( self.subreads, "all_subreads.fasta" )
        log.info('Extracting subreads from input files...')
        extract_subreads(args.input_file, 
                         self.subread_file,
                         min_length=args.min_read_length,
                         min_score=args.min_read_score,
                         max_count=args.max_count)
        check_output_file( self.subread_file )
        log.info('Finished the extraction of subread data\n')

    def align_new_subreads_to_contigs( self ):
        log.info('"Aligning all subreads to the HBAR contigs')
        self.contig_file = os.path.join( self.references, 'all_contigs.fasta' )
        self.subreads_to_contigs = os.path.join(self.alignments, 'subreads_to_contigs.m1')   
        subread_count = fasta_size( self.subread_file )
        contig_count = fasta_size( self.contig_file )
        # Run BLASR
        log.info("Aligning {0} contigs to the {1} contigs".format(subread_count, contig_count))
        blasr_args = {'nproc': args.nproc,
                      'out': self.subreads_to_contigs,
                      'bestn': 1,
                      'nCandidates': contig_count,
                      'noSplitSubreads': True}
        run_blasr( self.subread_file, 
                   self.contig_file,
                   blasr_args )
        # Check and save the output
        check_output_file( self.subreads_to_contigs )
        log.info("Finished aligning subreads to the created contigs\n")
        self.subread_contig_dict = create_m1_reference( self.subreads_to_contigs )

    def read_on_target_contigs( self ):
        log.info('Reading list of On-Target Contigs')
        self.on_target_contigs = os.path.join( self.references, 'hla_contig_ids.txt' )
        on_target_contig_ids = read_list_file( self.on_target_contigs )
        log.info('Finished reading the list of On-Target contigs\n')

    def separate_off_target_subreads( self ):
        log.info('Separating off-target subreads by contig alignment')
        self.hla_subreads = os.path.join(self.subreads, 'hla_subreads.fasta')
        self.non_hla_subreads = os.path.join(self.subreads, 'non_hla_subreads.fasta')
        separated_seqs = separate_aligned_sequences( self.subread_file, 
                                                     self.subread_contig_dict,
                                                     self.on_target_contigs )
        write_group(separated_seqs, 'selected', self.hla_subreads)
        write_group(separated_seqs, 'not_selected', self.non_hla_subreads)
        log.info('Finished separating off-target subreads\n')

    def create_ref_locus_dict( self ):    
        log.info('Creating a reference locus dict')
        self.reference_locus_dict = create_fofn_reference( args.reference_file )
        log.info('Finished creating reference locus dict\n')

    def create_phased_locus_dict( self ):    
        log.info('Creating phased contig locus dict')
        self.phased_to_reference = os.path.join(self.alignments, 'phased_to_reference.m5')   
        self.phased_reference_dict = create_m5_reference( self.phased_to_reference )
        self.phased_locus_dict = cross_ref_dict( self.phased_reference_dict, self.reference_locus_dict )
        log.info('Finished creating phased contig locus dict\n')

    def align_new_subreads_to_ref( self ):
        raw_output = os.path.join( self.results, 'Clusense', 'Final_Raw.fasta')
        reseq_output = os.path.join( self.results, 'Clusense', 'Final_Resequenced.fasta')
        if valid_file( reseq_output ):
            ref_file = reseq_output
        elif valid_file( raw_output ):
            ref_file = raw_output
        else:
            msg = 'No valid output reference file found!'
            log.error( msg )
            raise IOError( msg )
        log.info("Re-aligning HLA subreads to selected references")
        self.new_to_ref = os.path.join( self.alignments, "hla_subreads_to_ref.sam" )
        query_count = fasta_size( self.hla_subreads )
        ref_count = fasta_size( ref_file )
        log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
        blasr_args = {'nproc': args.nproc,
                      'noSplitSubreads': True,
                      'out': self.new_to_ref,
                      'sam': True,
                      'bestn': 1,
                      'nCandidates': ref_count}
        run_blasr( self.hla_subreads, 
                   ref_file, 
                   blasr_args )
        check_output_file( self.new_to_ref )
        self.subread_allele_dict = create_sam_reference( self.new_to_ref )
        self.subread_locus_dict = cross_ref_dict( self.subread_allele_dict, self.phased_locus_dict )
        log.info("Finished realigning HLA subreads to HLA contigs\n")

    def combine_subreads_by_allele(self):
        log.info("Separating subreads by their best matching allele")
        allele_fofn = os.path.join( self.subreads, "Alleles_Files.txt" )
        separated_seqs = separate_sequences( self.hla_subreads, 
                                             self.subread_allele_dict )
        allele_prefix = os.path.join(self.subreads, 'Allele')
        self.allele_files = write_all_groups( separated_seqs, allele_prefix )
        self.allele_files = [fn for fn in self.allele_files if not fn.endswith('Unmapped.fasta')]
        write_list_file( self.allele_files, allele_fofn )
        log.info('Finished separating subreads by Allele\n')

    def combine_subreads_by_locus(self):
        log.info("Separating subreads by the locus of their assigned contig")
        locus_fofn = os.path.join( self.subreads, "Locus_Files.txt" )
        separated_seqs = separate_sequences( self.hla_subreads, 
                                             self.subread_locus_dict )
        locus_prefix = os.path.join(self.subreads, 'Locus')
        self.locus_files = write_all_groups( separated_seqs, locus_prefix )
        self.locus_files = [fn for fn in self.locus_files if not fn.endswith('Unmapped.fasta')]
        write_list_file( self.locus_files, locus_fofn )
        log.info('Finished separating subreads by Locus\n')

if __name__ == '__main__':
    NewSubreadPipeline().run()
