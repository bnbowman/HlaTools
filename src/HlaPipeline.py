#! /usr/bin/env python
import os
import logging
import ConfigParser

from pbcore.io.FastaIO import FastaReader
from pbphase.clusense import Clusense
from pbphase.AmpliconAnalyzer import AmpliconAnalyzer

from pbhla.log import initialize_logger
from pbhla.fasta.update_orientation import update_orientation
from pbhla.amplicons.separate import separate_amplicons
from pbhla.arguments import args, parse_args
from pbhla.fofn import (create_baxh5_fofn,
                        write_sequence_fofn)
from pbhla.separate_sequences import (separate_sequences,
                                      separate_listed_sequences,
                                      separate_aligned_sequences)
from pbhla.fasta.utils import ( write_fasta,
                                fasta_size,
                                combine_fasta )
from pbhla.fasta.subset_sequences import subset_sequences
from pbhla.fasta.rename import rename_fofn
from pbhla.io.extract_subreads import extract_subreads
from pbhla.io.BlasrIO import add_header_to_m5
from pbhla.io.FofnIO import FofnReader
from pbhla.dictionary import ( create_m1_reference,
                               create_m5_reference,
                               create_amp_assem_reference,
                               create_sam_reference,
                               create_phased_reference,
                               filter_m5_file )
from pbhla.references.fofn import parse_reference_fofn
from pbhla.phasing.summarize import ( combine_clusense_output,
                                      summarize_amp_assem_output )
from pbhla.resequencing.identify import identify_resequencing_data
from pbhla.resequencing.SummarizeResequencing import combine_resequencing_output
from pbhla.annotation.summarize import ( summarize_contigs,
                                         summarize_resequenced,
                                         summarize_typing )
from pbhla.annotation.meta_summarize import meta_summarize_contigs
from pbhla.external.HbarTools import HbarRunner
from pbhla.external.commandline_tools import run_blasr
from pbhla.external.ClusterResequencer import ClusterResequencer
from pbhla.external.utils import align_best_reference
from pbhla.utils import *

log = logging.getLogger()

class HlaPipeline( object ):

    def __init__( self ):
        # Parse the options
        parse_args()
        self._config = ConfigParser.SafeConfigParser()
        self._config.read( args.config_file )


        # Initialize output folder and sub-folders
        self.subfolders = _initialize_folders( args.output )

        # Initialize logging
        log_file = os.path.join( args.output, "HLA_Pipeline.log" )
        initialize_logger( log_file )

    def __getattr__(self, item):
        if item in ['min_read_length', 'max_count']:
            return self._config.getint('DEFAULT', item)
        elif item in ['min_read_score', 'min_snr']:
            return self._config.getfloat('DEFAULT', item)
        elif item in ['clustering']:
            return self._config.getboolean('DEFAULT', item)
        else:
            return self._config.get('Global', item)

    def config(self, domain, item):
        hla_domain = 'HLA-%s' % domain
        print domain
        if domain in self._config.sections():
            print self._config.get( domain, item )
            return self._config.get( domain, item )
        elif hla_domain in self._config.sections():
            print self._config.get( hla_domain, item )
            return self._config.get( hla_domain, item )
        else:
            msg = 'Unrecognized configuration domain (%s)' % domain
            log.error( msg )
            raise KeyError( msg )

    def has_config(self, domain):
        hla_domain = 'HLA-%s' % domain
        if domain in self._config.sections():
            return True
        elif hla_domain in self._config.sections():
            return True
        else:
            return False

    def get_filepath(self, folder, filename):
        return os.path.join( self.subfolders[folder], filename )

    def run( self ):
        baxh5 = _create_baxh5_fofn( args.input_file, args.output )
        raw_subreads = self.extract_subread_data()

        # First we assemble the supplied data via HGAP / HBAR
        input_fofn = self.create_hbar_input_fofn( raw_subreads )
        contig_file = self.run_hbar_assembly( input_fofn )
        renamed_subreads = self.export_hbar_subreads()
        renaming_key = self.create_renaming_key( raw_subreads, renamed_subreads )

        # Second align the subreads and contigs to various references ...
        subread_contig_dict = self.align_subreads_to_contigs( renamed_subreads, contig_file )
        contig_genome_dict = self.align_contigs_to_genome( contig_file )
        hla_reference, metadata, loci = self.parse_reference()
        contig_reference_dict = self.align_contigs_to_reference( contig_file, hla_reference )
        contig_locus_dict = cross_ref_dict( contig_reference_dict, loci )
        subread_locus_dict = cross_ref_dict( subread_contig_dict, contig_locus_dict )

        # ... in order to separate the on-target from off-target sequences ...
        on_target = self.find_on_target_contigs( contig_genome_dict, contig_locus_dict )
        hla_contigs = self.separate_hla_contigs( contig_file, on_target )
        hla_subreads = self.separate_hla_subreads( renamed_subreads, subread_contig_dict, on_target )

        # ... and the different loci and contigs from each other
        contig_subread_fofn = self.separate_subreads_by_contig( hla_subreads, subread_contig_dict )
        locus_subread_fofn = self.separate_subreads_by_locus( hla_subreads, subread_locus_dict )

        # Rename and Phase each Locus
        renamed_subread_fofn = self.rename_subread_files( locus_subread_fofn, renaming_key )
        self.run_amp_analysis( renamed_subread_fofn )

        #self.rename_contig_subread_files()
        #self.phase_contigs_with_amp_assem()
        #self.summarize_amp_assem()
        #self.output_amp_assem_contigs()
        #self.align_amp_assem_to_reference()
        #self.align_subreads_to_amp_assem()
        #self.separate_subreads_by_allele()
        #self.summarize_amp_assem_results()
        #self.extract_final_amp_assem_contigs()
        #self.update_amp_assem_orientation()
        #self.align_final_amp_assem_to_reference()

    def extract_subread_data( self ):
        """
        Extract subreads meeting the DEFAULT requirements for analysis
        """
        log.info('Looking for raw subread data')
        # Dump all valid reads from the above files
        subread_file = self.get_filepath( "subreads", "all_subreads.fasta" )
        if valid_file( subread_file ):
            log.info("Using existing subread file\n")
            return subread_file

        log.info('No subread data found, extracting from input file(s)')
        extract_subreads( args.input_file, 
                          subread_file,
                          min_length=self.min_read_length,
                          min_score=self.min_read_score,
                          min_snr=self.min_snr,
                          max_count=self.max_count )
        check_output_file( subread_file )
        log.info('Finished extracting subread data from input\n')
        return subread_file

    def create_hbar_input_fofn( self, subread_file ):
        """
        Create a input FOFN for HBAR pointing to the raw subread data
        """
        log.info('Looking for HBAR input FOFN')
        input_fofn = self.get_filepath( 'HBAR', 'input.fofn' )
        if valid_file( input_fofn ):
            log.info('Using existing HBAR input FOFN file\n')
            return input_fofn

        log.info("No HBAR input FOFN found, creating from subread data")
        with open( input_fofn, 'w' ) as handle:
            handle.write( os.path.abspath( subread_file ))
        log.info('Finished creating input fofn\n')
        return input_fofn

    def run_hbar_assembly( self, input_fofn ):
        """
        Run HBAR to assemble rough contigs from the raw HLA subreads
        """
        log.info('Looking for HBAR-assembled HLA contigs')
        output = self.get_filepath( "HBAR", "3-CA/9-terminator/asm.utg.fasta" )
        contig_file = self.get_filepath( "references", "all_contigs.fasta" )
        if valid_file( output ):
            log.info("Using existing HBAR contig file")
        else: # Run HGAP
            log.info("No HBAR contig file found, initializing HbarRunner")
            hbar = HbarRunner( input_fofn,
                               self.subfolders["HBAR"],
                               min_length=self.min_read_length,
                               min_score=self.min_read_score )
            hbar()

        # Copy the contig file to a more convenient location
        check_output_file( output )
        copy_file( output, contig_file )
        log.info('Finished assembling subreads data with HBAR\n')
        return contig_file

    def export_hbar_subreads(self):
        """
        Export the HBAR-renamed subread data for downstream analysis
        """
        log.info('Looking for exported renamed subreads from HBAR')
        renamed_subreads = self.get_filepath( 'subreads', 'renamed_subreads.fasta' )
        if valid_file( renamed_subreads ):
            log.info('Using existing HBAR-exported subread file\n')
            return renamed_subreads

        log.info("No renamed subread file found, exporting from HBAR output...")
        fasta_folder = self.get_filepath( 'HBAR', '0-fasta_files' )
        for entry in os.listdir( fasta_folder ):
            if entry.endswith('_q.fa'):
                hbar_fasta = os.path.join( fasta_folder, entry )
                copy_file( hbar_fasta, renamed_subreads )
        check_output_file( renamed_subreads )
        log.info('Finished exporting HBAR-renamed subreads\n')
        return renamed_subreads

    def create_renaming_key(self, raw_subreads, renamed_subreads ):
        """
        Create a key for translating HBAR subread names to canonical PacBio names
        """
        log.info("Looking for Raw<--->HBAR subread renaming key")
        renaming_key = self.get_filepath( 'subreads', 'renaming_key.txt' )
        if valid_file( renaming_key ):
            log.info('Using existing subread renaming key\n')
            return renaming_key

        log.info("No subread renaming key round, creating one...")
        # Compare the two files to make sure they're equivalent
        raw_count = fasta_size( raw_subreads )
        new_count = fasta_size( renamed_subreads )
        try:
            assert raw_count == new_count 
        except AssertionError:
            msg = 'The number of raw subreads (%s) does not ' % raw_count + \
                  'match the number of renamed reads (%s)' % new_count
            log.info( msg )
            raise ValueError( msg )
        # Write out the pairs of names to file
        with open( renaming_key, 'w') as handle:
            for raw, renamed in zip( FastaReader(raw_subreads), FastaReader(renamed_subreads) ):
                raw_name = raw.name.split()[0]
                new_name = renamed.name.split()[0]
                handle.write('%s\t%s\n' % (new_name, raw_name))
        check_output_file( renaming_key )
        log.info("Finished creating subread renaming key\n")
        return renaming_key

    def align_subreads_to_contigs(self, subread_file, contig_file ):
        """
        Align the subreads to the contigs assembled by HBAR
        """
        log.info("Looking for Subread-to-Contig alignment data")
        subread_contig_align = self.get_filepath( 'alignments', 'subreads_to_contigs.m1' )
        if valid_file( subread_contig_align ):
            log.info("Using existing Subread->Contig alignment file\n")
        else:
            log.info("No Subread->Contig alignment found, creating...")
            align_best_reference( subread_file, contig_file, output=subread_contig_align )
            check_output_file( subread_contig_align )
            log.info("Finished aligning subreads to the HBAR contigs\n")
        return create_m1_reference( subread_contig_align )

    def align_contigs_to_genome(self, contig_file):
        log.info("Looking for Contig-to-Genome alignment data")
        contig_genome_align = self.get_filepath( 'alignments', 'contigs_to_genome.m1' )
        if valid_file( contig_genome_align ):
            log.info("Using existing Contig->Genome alignment file\n")
        else:
            log.info("No Contig->Genome alignment found, creating...")
            align_best_reference( contig_file, self.human_reference, output=contig_genome_align )
            check_output_file( contig_genome_align )
            log.info("Finished aligning contigs to the genomic reference\n")
        return create_m1_reference( contig_genome_align )

    def parse_reference(self):
        """
        Parse HLA data from the configured reference FOFN
        """
        log.info("Parsing the supplied FOFN of HLA reference data")
        hla_reference_seqs = self.get_filepath( "references", "HLA_references.fasta" )
        sequences, metadata, loci = parse_reference_fofn( self.hla_reference )
        log.info("Writing collected HLA reference sequences to file")
        write_fasta( sequences, hla_reference_seqs )
        check_output_file( hla_reference_seqs )
        log.info("Finished parsing the HLA reference data\n")
        return hla_reference_seqs, metadata, loci
    
    def align_contigs_to_reference(self, contig_file, reference_file):
        """
        Align HBAR contigs to an HLA reference Fasta
        """
        log.info("Looking for Contig-to-Reference alignment data")
        contig_reference_align = self.get_filepath( 'alignments', 'contigs_to_reference.m1' )
        if valid_file( contig_reference_align ):
            log.info("Using an existing Contig->Reference alignment file\n")
        else:
            log.info("No Contig->Reference alignment found, creating...")
            align_best_reference( contig_file, reference_file, output=contig_reference_align )
            check_output_file( contig_reference_align )
            log.info("Finished aligning contigs to the HLA reference data\n")
        return create_m1_reference( contig_reference_align )

    def find_on_target_contigs( self, contig_genome, contig_loci ):
        """
        Identify on-target contigs based on their genomic and reference alignments
        """
        log.info('Identifying on-target contigs from alignment data')
        on_target_contig_ids = [c for c in contig_loci if contig_genome[c] == 'chr6']
        log.info('Finished identifying on-target contigs\n')
        return on_target_contig_ids

    def separate_hla_contigs( self, contig_file, on_target_ids ):
        """
        Separate the ON/OFF target contig sequences from each other
        """
        log.info('Looking for separated on/off target contig data')
        hla_contigs = self.get_filepath( 'references', 'hla_contigs.fasta' )
        other_contigs = self.get_filepath( 'references', 'other_contigs.fasta' )
        if valid_file( hla_contigs ):
            log.info("Using existing separated contig files\n")
            return hla_contigs
        else:
            log.info('No separated contig files found, creating...')
            separate_listed_sequences( contig_file, on_target_ids, hla_contigs, other_contigs )
            log.info('Finished separating on/off-target contigs\n')
        return hla_contigs
    
    def separate_hla_subreads(self, subread_file, subread_contigs, on_target_ids ):
        """
        Separate the ON/OFF target subreads based on their aligned contig
        """
        log.info('Looking for separated on/off target subread data')
        hla_subreads = self.get_filepath( 'subreads', 'hla_subreads.fasta' )
        other_subreads = self.get_filepath( 'subreads', 'other_subreads.fasta' )
        if valid_file( hla_subreads ):
            log.info("Using existing separated subread files\n")
            return hla_subreads
        else:
            log.info('No separated subread files found, creating...')
            separate_aligned_sequences( subread_file, subread_contigs, on_target_ids,
                                        hla_subreads, other_subreads )
            log.info('Finished separating on/off-target subreads\n')
        return hla_subreads

    def separate_subreads_by_contig( self, hla_subreads, subread_contigs ):
        """
        Separate the on-target subreads by the HBAR contig to which they align
        """
        log.info("Looking for contig-specific subread files")
        contig_subread_fofn = self.get_filepath( "subreads", "Contig_Subread_Files.fofn" )
        if valid_file( contig_subread_fofn ):
            log.info("Using existing contig subread files and FOFN\n")
            return contig_subread_fofn
        log.info("No contig subread files found, creating...")
        subread_prefix = self.get_filepath( "subreads", "Contig" )
        contig_files = separate_sequences( hla_subreads, subread_contigs, subread_prefix )
        write_sequence_fofn( contig_files, contig_subread_fofn )
        log.info('Finished separating subreads by contig\n')
        return contig_subread_fofn

    def separate_subreads_by_locus( self, hla_subreads, subread_loci ):
        """
        Separate the on-target subreads by the HBAR contig to which they align
        """
        log.info("Looking for locus-specific subread files")
        locus_subread_fofn = self.get_filepath( "subreads", "Locus_Subread_Files.fofn" )
        if valid_file( locus_subread_fofn ):
            log.info("Using existing locus subread files and FOFN\n")
            return locus_subread_fofn
        log.info("No locus subread files found, creating...")
        subread_prefix = self.get_filepath( "subreads", "Locus" )
        locus_files = separate_sequences( hla_subreads, subread_loci, subread_prefix )
        write_sequence_fofn( locus_files, locus_subread_fofn )
        log.info('Finished separating subreads by locus\n')
        return locus_subread_fofn

    def rename_subread_files( self, subread_fofn, renaming_key ):
        log.info("Looking for FOFN of renamed subread files")
        renamed_fofn = '.'.join(subread_fofn.split('.')[:-1]) + '_renamed.fofn'
        if valid_file( renamed_fofn ):
            log.info("Using existing FOFN of renamed subread files\n")
            return renamed_fofn
        log.info("No renamed FOFN found, creating...")
        rename_fofn( subread_fofn, renamed_fofn, renaming_key )
        log.info("Finished renaming subread files\n")
        return renamed_fofn

    def run_amp_analysis( self, renamed_fofn ):
        log.info("Looking for phased AmpliconAnalysis results")
        analyzer = AmpliconAnalyzer( self.smrt_path, args.nproc )
        for subread_file in FofnReader( renamed_fofn ):

            # Check if the source of the current file has a configuration
            source = get_file_source( subread_file )
            output_folder = self.get_filepath( 'phasing', source )
            if self.has_config( source ):
                log.info('Found an AmpliconAnalysis configuration for "%s"' % source)
            else:
                log.info('No AmpliconAnalysis configuration found for "%s"' % source)
                continue

            # For sources with a configuration, run AA if output DNE
            if os.path.exists( output_folder ):
                log.info('AmpliconAnalyzer output detected for "%s", skipping...' % source)
            else:
                log.info('Phasing subreads from "%s"' % source)
                analyzer_args = {'whiteList': subread_file,
                                 'sampleName': '_' + source,
                                 'minLength': self.config(source, 'min_read_length'),
                                 'minReadScore': self.config(source, 'min_read_score'),
                                 'minSnr': self.config(source, 'min_snr'),
                                 'noClustering': self.config(source, 'disable_clustering')}
                print analyzer_args
                analyzer.run( args.input_file, output_folder, analyzer_args )
        log.info('Finished phasing subreads with Amplicon Analysis\n')

    def summarize_amp_assem( self ):
        log.info("Summarizing the output from Amplicon Assembly")
        output_folder = os.path.join(self.phasing_results, 'AmpliconAnalyzer')
        self.amp_assem_results = summarize_amp_assem_output( self.amp_assem, 
                                                             output_folder )
        log.info('Finished Phasing subreads with Amplicon Assembly')

    def output_amp_assem_contigs( self ):
        log.info("Combining the output from Amplicon Assembly")
        self.amp_assem_contigs = os.path.join( self.references, 'amp_assem_contigs.fasta')
        result_files = read_list_file( self.amp_assem_results )
        combine_fasta( result_files, self.amp_assem_contigs )
        log.info("Finished combining the output from Amplicon Assembly")

    def align_amp_assem_to_reference(self):
        log.info('Aligning all phased HLA contigs to the Reference sequences')
        self.amp_assem_to_reference = os.path.join(self.alignments, 'amp_assem_to_reference.m1')   
        contig_count = fasta_size( self.amp_assem_contigs )
        reference_count = fasta_size( self.reference_seqs )
        log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
        # Run BLASR
        blasr_args = {'nproc': args.nproc,
                      'noSplitSubreads': True,
                      'out': self.amp_assem_to_reference,
                      'm': 1,
                      'bestn': 1,
                      'nCandidates': reference_count}
        run_blasr( self.amp_assem_contigs, 
                   self.reference_seqs, 
                   blasr_args )
        check_output_file( self.amp_assem_to_reference )
        # Check and save the output
        self.amp_assem_reference_dict = create_m1_reference( self.amp_assem_to_reference )
        self.amp_assem_locus_dict = create_amp_assem_reference( self.amp_assem_to_reference )
        #self.subread_locus_dict = cross_ref_dict( self.phased_read_dict, self.phased_locus_dict )
        log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def align_subreads_to_amp_assem(self):
        log.info("Re-aligning HLA subreads to selected references")
        self.hla_to_amp_assem = os.path.join( self.alignments, "hla_subreads_to_amp_assem.sam" )
        if valid_file( self.hla_to_amp_assem ):
            log.info('Found existing SAM file "%s"' % self.hla_to_amp_assem)
            log.info("Skipping realignment step...\n")
        else:
            query_count = fasta_size( self.hla_subreads )
            ref_count = fasta_size( self.amp_assem_contigs )
            log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
            blasr_args = {'nproc': args.nproc,
                          'noSplitSubreads': True,
                          'out': self.hla_to_amp_assem,
                          'sam': True,
                          'bestn': 1,
                          'nCandidates': ref_count}
            run_blasr( self.hla_subreads, 
                       self.amp_assem_contigs, 
                       blasr_args )
            check_output_file( self.hla_to_amp_assem )
        self.subread_allele_dict = create_sam_reference( self.hla_to_amp_assem )
        self.subread_locus_dict = cross_ref_dict( self.subread_allele_dict, 
                                                  self.amp_assem_locus_dict )
        log.info("Finished realigning HLA subreads to HLA contigs\n")

    def separate_subreads_by_allele(self):
        log.info("Separating subreads by their best matching allele")
        self.allele_subread_fofn = os.path.join( self.subreads, "Alleles_Subread_Files.fofn" )
        if valid_file( self.allele_subread_fofn ):
            log.info('Found existing Allele File List "%s"' % self.allele_subread_fofn)
            log.info("Skipping subread separation step...\n")
            self.allele_files = read_list_file( self.allele_subread_fofn )
            return
        separated_seqs = separate_sequences( self.hla_subreads, 
                                             self.subread_allele_dict )
        allele_prefix = os.path.join(self.subreads, 'Allele')
        self.allele_files = write_all_groups( separated_seqs, allele_prefix )
        self.allele_files = [fn for fn in self.allele_files if not fn.endswith('Unmapped.fasta')]
        write_list_file( self.allele_files, self.allele_subread_fofn )
        log.info('Finished separating subreads by Allele')

    def summarize_amp_assem_results(self):
        log.info("Picking the best contigs from the phasing results")
        self.amp_assem_results = os.path.join( self.results, 'AmpliconAssembly' )
        create_directory( self.amp_assem_results )
        log.info("Summarizing individual loci")
        summaries = summarize_contigs( self.amp_assem_contigs, 
                                       self.allele_subread_fofn, 
                                       self.amp_assem_locus_dict, 
                                       self.amp_assem_to_reference,
                                       self.amp_assem_results )
        log.info("Combining the summaries of the various HLA loci")
        self.meta_summary = os.path.join( self.amp_assem_results, 'Locus_Calls.txt')
        meta_summarize_contigs( summaries, 
                                self.reference_metadata,
                                self.meta_summary,
                                excluded=args.exclude)
        log.info("Finished selected best contigs")

    def extract_final_amp_assem_contigs(self):
        log.info('Extracting the best contigs to their own file')
        self.final_contigs = os.path.join(self.references, 'final_contigs.fasta')
        subset_sequences( self.amp_assem_contigs,
                          self.meta_summary,
                          self.final_contigs )
        log.info('Finished extracting the best resequenced contigs\n')

    def update_amp_assem_orientation(self):
        log.info('Converting all resequenced contigs to the forward orientation')
        self.reoriented = os.path.join(self.references, 'reoriented_final_contigs.fasta')
        self.final_output = os.path.join(self.amp_assem_results, 'Final_Sequences.fasta')
        update_orientation( self.final_contigs, 
                            self.amp_assem_to_reference,
                            self.reoriented )
        shutil.copy( self.reoriented, self.final_output )
        log.info('Finished updating the orientation of the resequenced contigs')

    def align_final_amp_assem_to_reference(self):
        log.info('Aligning all final phased HLA contigs to the reference sequences')
        self.final_to_reference = os.path.join(self.amp_assem_results, 'Final_to_Reference.m5')   
        contig_count = fasta_size( self.amp_assem_contigs )
        reference_count = fasta_size( self.reference_seqs )
        log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
        # Run BLASR
        blasr_args = {'nproc': args.nproc,
                      'noSplitSubreads': True,
                      'out': self.final_to_reference,
                      'm': 5,
                      'bestn': 1,
                      'nCandidates': reference_count}
        run_blasr( self.final_output, 
                   self.reference_seqs, 
                   blasr_args )
        check_output_file( self.final_to_reference )
        add_header_to_m5( self.final_to_reference )
        check_output_file( self.final_to_reference )
        log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def separate_contigs(self):
        log.info("Separating remaining contigs into individual files")
        contig_fofn = os.path.join( self.references, "Contig_Files.txt" )
        if valid_file( contig_fofn ):
            log.info('Found existing Contig File List "{0}"'.format(contig_fofn))
            log.info("Skipping contig separation step...\n")
            self.contig_files = read_list_file( contig_fofn )
            return
        separated_seqs = separate_sequences( self.selected_contigs )
        contig_prefix = os.path.join(self.references, 'Contig')
        self.contig_files = write_all_groups( separated_seqs, contig_prefix )
        write_list_file( self.contig_files, contig_fofn)
        log.info('Finished separating contigs\n')

    def phase_reads_with_clusense(self):
        log.info("Phasing subreads with Clusense")
        self.clusense = os.path.join( self.phasing, 'Clusense' )
        create_directory( self.clusense )
        for subread_fn, contig_fn in zip(self.subread_files, self.contig_files):
            folder_name = os.path.basename(contig_fn).split('.')[0]
            output_folder = os.path.join(self.clusense, folder_name)
            if os.path.exists( output_folder ):
                log.info('"Clusense output detected for "{0}", skipping...'.format(folder_name))
            else:
                log.info("Phasing subreads for {0}".format(folder_name))
                Clusense(subread_fn, contig_fn, output_folder, threshold=0.1)
            cleanup_directory( output_folder )
        log.info('Finished Phasing subreads with Clusense')

    def summarize_clusense(self):
        log.info("Summarizing the output from Clusense")
        output_folder = os.path.join(self.phasing_results, 'Clusense')
        self.clusense_cns_fofn, self.clusense_read_fofn = combine_clusense_output(self.clusense, output_folder)
        self.phased_read_dict = create_phased_reference( self.clusense_read_fofn )
        log.info('Finished Phasing subreads with Clusense')

    def output_phased_contigs(self):
        log.info("Creating a combined reference of non-HLA and phased-HLA contigs")
        self.phased_contigs = os.path.join(self.references, 'phased_contigs.fasta')
        clusense_cns_files = read_list_file( self.clusense_cns_fofn )
        fasta_files = list(clusense_cns_files) + [self.non_hla_contigs]
        combine_fasta(fasta_files, self.phased_contigs)
        self.phased_hla = os.path.join(self.references, 'phased_hla.fasta')
        combine_fasta(clusense_cns_files, self.phased_hla)
        log.info('Finished creating combined reference')

    def align_phased_to_reference(self):
        log.info('Aligning all phased HLA contigs to the Reference sequences')
        raw_phased_to_reference = os.path.join(self.alignments, 'raw_phased_to_reference.m5')   
        self.phased_to_reference = os.path.join(self.alignments, 'phased_to_reference.m5')   
        if valid_file( self.phased_to_reference ):
            log.info('Found existing alignment file "{0}"'.format(self.phased_to_reference))
            log.info("Skipping alignment step...\n")
        else:
            contig_count = fasta_size( self.phased_hla )
            reference_count = fasta_size( self.reference_seqs )
            # Run BLASR
            log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
            blasr_args = {'nproc': args.nproc,
                          'noSplitSubreads': True,
                          'out': raw_phased_to_reference,
                          'm': 5,
                          'bestn': 5,
                          'nCandidates': reference_count}
            run_blasr( self.phased_hla, 
                       self.reference_seqs, 
                       blasr_args )
            filter_m5_file( raw_phased_to_reference, self.phased_to_reference )
            # Check and save the output
            check_output_file( self.phased_to_reference )
        self.phased_reference_dict = create_m5_reference( self.phased_to_reference )
        self.phased_locus_dict = cross_ref_dict( self.phased_reference_dict, self.reference_locus_dict )
        self.subread_locus_dict = cross_ref_dict( self.phased_read_dict, self.phased_locus_dict )
        log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def combine_subreads_by_locus(self):
        log.info("Separating subreads by the locus of their assigned contig")
        locus_fofn = os.path.join( self.subreads, "Locus_Files.txt" )
        if valid_file( locus_fofn ):
            log.info('Found existing Locus File List "{0}"'.format(locus_fofn))
            log.info("Skipping subread separation step...\n")
            self.locus_files = read_list_file( locus_fofn )
            return
        separated_seqs = separate_sequences( self.hla_subreads, 
                                             self.subread_locus_dict )
        locus_prefix = os.path.join(self.subreads, 'Locus')
        self.locus_files = write_all_groups( separated_seqs, locus_prefix )
        self.locus_files = [fn for fn in self.locus_files if not fn.endswith('Unmapped.fasta')]
        write_list_file( self.locus_files, locus_fofn )
        log.info('Finished separating subreads by Locus')

    def summarize_phased_by_locus(self):
        log.info("Picking the best contigs from the phasing results")
        self.clusense_results = os.path.join( self.results, 'Clusense' )
        create_directory( self.clusense_results )
        log.info("Summarizing individual loci")
        summaries = summarize_contigs( self.phased_hla, 
                                       self.clusense_read_fofn, 
                                       self.phased_locus_dict, 
                                       self.phased_to_reference,
                                       self.clusense_results )
        log.info("Combining the summaries of the various HLA loci")
        self.meta_summary = os.path.join( self.clusense_results, 'Locus_Calls.txt')
        meta_summarize_contigs( summaries, 
                                self.meta_summary,
                                excluded=args.exclude)
        log.info("Finished selected best contigs")

    def extract_best_contigs(self):
        log.info('Extracting the best contigs to their own file')
        self.raw_output = os.path.join(self.clusense_results, 'Final_Raw.fasta')
        if valid_file( self.raw_output ):
            log.info('Found existing Resequencing Output')
            log.info('Skipping...\n')
            return
        subset_sequences( self.phased_hla,
                          self.meta_summary,
                          self.raw_output )
        log.info('Finished extracting the best resequenced contigs\n')

    def identify_resequencing_files(self):
        log.info("Identifying files for resequencing based on selected contigs")
        self.resequencing_pairs = identify_resequencing_data( self.meta_summary,
                                                              self.clusense_cns_fofn, 
                                                              self.clusense_read_fofn )
        log.info('Finished the idenifying the requisite files')

    def resequence_contigs(self):
        log.info('Resequencing selected contigs')
        self.resequencing_dir = os.path.join( self.resequencing, 'Clusense' )
        create_directory( self.resequencing_dir )
        for cluster_name, cns_file, read_file in self.resequencing_pairs:
            output_folder = os.path.join( self.resequencing_dir, cluster_name )
            if os.path.exists( output_folder ):
                log.info('Resequencing output detected for "{0}", skipping...'.format(cluster_name))
            else:
                log.info('Resequencing cluster "{0}"'.format(cluster_name))
                resequencer = ClusterResequencer(read_file,
                                                 cns_file,
                                                 self.baxh5_fofn,
                                                 setup=args.smrt_path,
                                                 names=self.renaming_key,
                                                 output=output_folder,
                                                 nproc=args.nproc)
                resequencer()
        log.info('Finished resequencing selected contigs')

    def summarize_resequencing(self):
        log.info("Summarizing the output from Resequencing")
        output_folder = os.path.join(self.resequencing_results, 'Clusense')
        self.resequencing_fasta_fofn, self.resequencing_fastq_fofn = combine_resequencing_output(self.resequencing_dir, output_folder)
        log.info('Finished summarizing the resequencing results')

    def output_resequenced_contigs(self):
        log.info("Creating a combined reference of the selected and resequenced HLA contigs")
        self.resequenced_hla_contigs = os.path.join(self.references, 'resequenced_hla_contigs.fasta')
        fasta_files = read_list_file( self.resequencing_fasta_fofn )
        combine_fasta(fasta_files, self.resequenced_hla_contigs)
        log.info('Finished creating combined reference')

    def align_resequenced_to_reference(self):
        log.info('"Aligning all resequenced HLA contigs to the Reference sequences')
        self.resequenced_to_reference = os.path.join(self.alignments, 'resequenced_to_reference.m1')   
        if valid_file( self.resequenced_to_reference ):
            log.info('Found existing alignment file "{0}"'.format(self.resequenced_to_reference))
            log.info("Skipping alignment step...\n")
        else:
            contig_count = fasta_size( self.resequenced_hla_contigs )
            reference_count = fasta_size( self.reference_seqs )
            # Run BLASR
            log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
            blasr_args = {'nproc': args.nproc,
                          'noSplitSubreads': True,
                          'out': self.resequenced_to_reference,
                          'bestn': 1,
                          'nCandidates': reference_count}
            run_blasr( self.resequenced_hla_contigs, 
                       self.reference_seqs, 
                       blasr_args )
            # Check and save the output
            check_output_file( self.resequenced_to_reference )
        log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def add_resequencing_summary(self):
        log.info('Summarizing the the resequenced contigs')
        self.resequenced_summary = os.path.join(self.clusense_results, 'Resequenced_Calls.txt')
        if valid_file( self.resequenced_summary ):
            log.info('Found existing resequenced summary file "%s"' % self.resequenced_summary )
            log.info('Skipping...\n')
            return
        summarize_resequenced( self.meta_summary, 
                               self.resequenced_to_reference, 
                               self.resequenced_summary )
        log.info('Finished summarizing the resequenced contigs\n')

    def extract_best_resequenced_contigs(self):
        log.info('Extracting the best resequenced contigs to their own file')
        self.resequenced_output = os.path.join(self.clusense_results, 'Final_Resequenced.fasta')
        if valid_file( self.resequenced_output ):
            log.info('Found existing Resequencing Output')
            log.info('Skipping...\n')
            return
        subset_sequences( self.resequenced_hla_contigs,
                          self.resequenced_summary,
                          self.resequenced_output )
        log.info('Finished extracting the best resequenced contigs\n')

    def align_subreads_to_raw(self):
        log.info("Re-aligning HLA subreads to selected references")
        self.hla_to_raw = os.path.join( self.alignments, "hla_subreads_to_raw.sam" )
        if valid_file( self.hla_to_raw ):
            log.info('Found existing SAM file "%s"' % self.hla_to_raw)
            log.info("Skipping realignment step...\n")
        else:
            query_count = fasta_size( self.hla_subreads )
            ref_count = fasta_size( self.raw_output )
            log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
            blasr_args = {'nproc': args.nproc,
                          'noSplitSubreads': True,
                          'out': self.hla_to_raw,
                          'sam': True,
                          'bestn': 1,
                          'nCandidates': ref_count}
            run_blasr( self.hla_subreads, 
                       self.raw_output, 
                       blasr_args )
            check_output_file( self.hla_to_raw )
        self.subread_allele_dict = create_sam_reference( self.hla_to_raw )
        self.subread_locus_dict = cross_ref_dict( self.subread_allele_dict, self.phased_locus_dict )
        log.info("Finished realigning HLA subreads to HLA contigs\n")

    def align_subreads_to_resequenced(self):
        log.info("Re-aligning HLA subreads to selected references")
        self.hla_to_resequenced = os.path.join( self.alignments, "hla_subreads_to_resequenced.sam" )
        if valid_file( self.hla_to_resequenced ):
            log.info('Found existing SAM file "%s"' % self.hla_to_resequenced)
            log.info("Skipping realignment step...\n")
        else:
            query_count = fasta_size( self.hla_subreads )
            ref_count = fasta_size( self.resequenced_output )
            log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
            blasr_args = {'nproc': args.nproc,
                          'noSplitSubreads': True,
                          'out': self.hla_to_resequenced,
                          'sam': True,
                          'bestn': 1,
                          'nCandidates': ref_count}
            run_blasr( self.hla_subreads, 
                       self.resequenced_output, 
                       blasr_args)
            check_output_file( self.hla_to_resequenced )
        self.subread_allele_dict = create_sam_reference( self.hla_to_resequenced )
        self.subread_locus_dict = cross_ref_dict( self.subread_allele_dict, self.phased_locus_dict )
        log.info("Finished realigning HLA subreads to HLA contigs\n")

    def combine_subreads_by_allele(self):
        log.info("Separating subreads by their best matching allele")
        allele_fofn = os.path.join( self.subreads, "Alleles_Files.txt" )
        if valid_file( allele_fofn ):
            log.info('Found existing Allele File List "%s"' % allele_fofn)
            log.info("Skipping subread separation step...\n")
            self.allele_files = read_list_file( allele_fofn )
            return
        separated_seqs = separate_sequences( self.hla_subreads, 
                                             self.subread_allele_dict )
        allele_prefix = os.path.join(self.subreads, 'Allele')
        self.allele_files = write_all_groups( separated_seqs, allele_prefix )
        self.allele_files = [fn for fn in self.allele_files if not fn.endswith('Unmapped.fasta')]
        write_list_file( self.allele_files, allele_fofn )
        log.info('Finished separating subreads by Allele')

    def combine_subreads_by_locus(self):
        log.info("Separating subreads by the locus of their assigned contig")
        locus_fofn = os.path.join( self.subreads, "Locus_Files.txt" )
        if valid_file( locus_fofn ):
            log.info('Found existing Locus File List "{0}"'.format(locus_fofn))
            log.info("Skipping subread separation step...\n")
            self.locus_files = read_list_file( locus_fofn )
            return
        separated_seqs = separate_sequences( self.hla_subreads, 
                                             self.subread_locus_dict )
        locus_prefix = os.path.join(self.subreads, 'Locus')
        self.locus_files = write_all_groups( separated_seqs, locus_prefix )
        self.locus_files = [fn for fn in self.locus_files if not fn.endswith('Unmapped.fasta')]
        write_list_file( self.locus_files, locus_fofn )
        log.info('Finished separating subreads by Locus')

    def type_hla_sequences(self):
        log.info('Typing the selected HLA consensus sequences')
        if args.msa is None:
            log.info("No HLA MSA information supplied, can't perform HLA typing")
            return 
        self.gdna_types = os.path.join( self.annotation, 'gDNA_allele_calls.txt' )
        self.cdna_types = os.path.join( self.annotation, 'cDNA_allele_calls.txt' )
        type_hla( args.msa, 
                  self.reoriented, 
                  self.phased_locus_dict,
                  self.annotation )
        log.info('Finished typing the selected HLA sequences\n')

    def summarize_hla_typings(self):
        log.info('Typing the selected HLA consensus sequences')
        self.typing_summary = os.path.join(self.clusense_results, 'Final_Summary.txt')
        summarize_typing( self.resequenced_summary,
                          self.gdna_types,
                          self.cdna_types,
                          self.typing_summary )
        log.info('Finished typing the selected HLA sequences\n')


def _initialize_folders( output ):
    """
    Create the Main and Sub Output directories
    """
    # Create Main
    log.info("Creating output directories")
    output = os.path.abspath( output )
    create_directory( output )

    # Create sub-directories
    subfolders = {}
    for dir_name in ['HBAR', 'references', 'subreads', 'results',
                     'alignments', 'phasing', 'phasing_results']:
        sub_dir = os.path.join( output, dir_name )
        create_directory( sub_dir )
        subfolders[dir_name] = sub_dir
    return subfolders

def _create_baxh5_fofn( input_file, output ):
    """
    Convert any BasH5 input files to BaxH5 to avoid file-type problems
    """
    log.info('Creating FOFN of Bax.H5 files')
    baxh5_fofn = os.path.join( output, 'baxh5.fofn' )
    if valid_file( baxh5_fofn ):
        log.info("Using existing Bax.H5 FOFN file")
        return baxh5_fofn

    log.info("No existing Bax.H5 fofn found")
    create_baxh5_fofn( input_file, baxh5_fofn )
    check_output_file( baxh5_fofn )
    log.info('Finished writing Bax.H5 fofn file\n')
    return baxh5_fofn

if __name__ == '__main__':
    HlaPipeline().run()
