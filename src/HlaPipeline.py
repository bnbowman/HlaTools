#! /usr/bin/env python
import os
import logging
import ConfigParser

from pbcore.io.FastaIO import FastaReader
from pbphase.AmpliconAnalyzer import AmpliconAnalyzer

from pbhla.log import initialize_logger
from pbhla.arguments import args, parse_args
from pbhla.fofn import ( create_baxh5_fofn,
                         write_sequence_fofn )
from pbhla.separate_sequences import ( separate_sequences,
                                       separate_listed_sequences,
                                       separate_aligned_sequences )
from pbhla.fasta.utils import ( write_fasta,
                                fasta_size )
from pbhla.fasta.rename import rename_fofn
from pbhla.io.extract_subreads import extract_subreads
from pbhla.io.FofnIO import FofnReader
from pbhla.amplicon_analysis.chimeras import ChimeraDetector
from pbhla.dictionary import create_m1_reference
from pbhla.references.fofn import parse_reference_fofn
from pbhla.phasing.combine import ( combine_amp_analysis,
                                    combine_resequencing,
                                    combine_fastq )
from pbhla.external.HbarTools import HbarRunner
from pbhla.external.utils import align_best_reference
from pbhla.resequencing.Resequencer import Resequencer
from pbhla.utilities.rename_fastq import rename_resequencing
from pbhla.utilities.filter_fastq import filter_fastq
from pbhla.typing.sequences import type_sequences
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
        initialize_logger( log, log_file=log_file )

    def __getattr__(self, item):
        if item in ['min_read_length', 'max_read_length', 'min_num_reads', 'min_consensus_length', 'max_count']:
            return self._config.getint('DEFAULT', item)
        elif item in ['min_read_score', 'min_snr']:
            return self._config.getfloat('DEFAULT', item)
        elif item in ['clustering']:
            return self._config.getboolean('DEFAULT', item)
        else:
            return self._config.get('Global', item)

    def config(self, domain, item):
        domain = self._check_domain( domain )
        print domain, item
        return self._config.get( domain, item )

    def _check_domain(self, domain):
        hla_domain = 'HLA-%s' % domain
        numbered_domain = hla_domain + '1'
        if domain in self._config.sections():
            return domain
        elif hla_domain in self._config.sections():
            return hla_domain
        elif numbered_domain in self._config.sections():
            return numbered_domain
        else:
            return None

    def to_be_phased(self, domain):
        domain = self._check_domain( domain )
        if domain:
            if self._config.getboolean( domain, 'use_amp_analysis' ):
                log.info("AmpliconAnalysis enabled for %s" % domain)
                return True
            else:
                log.info("AmpliconAnalysis disabled for %s" % domain)
                return False
        log.info("No configuration detected for %s" % domain)
        return False

    def to_be_resequenced(self, domain):
        domain = self._check_domain( domain )
        print "Reseq", domain
        if domain:
            if self._config.getboolean( domain, 'use_resequencing' ):
                log.info("Resequencing enabled for %s" % domain)
                return True
            else:
                log.info("resequencing disabled for %s" % domain)
                return False
        log.info("No configuration detected for %s" % domain)
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
        locus_contig_fofn = self.separate_contigs_by_locus( hla_contigs, contig_locus_dict )

        # Rename the subreads for each locus
        renamed_subread_fofn = self.rename_subread_files( locus_subread_fofn, renaming_key )

        # Use the renamed subreads to Phase with AA and remove AA-chimeras
        self.run_amp_analysis( renamed_subread_fofn )
        phasing_results = self.combine_phasing_results()
        #good_results = self.remove_chimeras( phasing_results )

        # Use the renamed subreads to resequence any other loci/HBAR contigs
        self.run_resequencing( baxh5, renamed_subread_fofn, locus_contig_fofn )
        reseq_results = self.combine_resequencing_results()
        renamed_results = self.rename_resequencing_results( reseq_results )

        # Combine the results from AA and resequencing and type the results
        combined_results = self.combine_phasing_and_reseq( phasing_results, renamed_results )
        filtered_results = self.filter_combined_results( combined_results )
        typing_sequences = self.copy_sequences_for_typing( filtered_results )
        self.type_hla_sequences( typing_sequences )

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
                          max_length=self.max_read_length,
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
        Separate the on-target subreads by the locus to which their HBAR contig aligns
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

    def separate_contigs_by_locus( self, hla_contigs, contig_loci ):
        """
        Separate the on-target contigs by the locus to which they align
        """
        log.info("Looking for locus-specific contig files")
        locus_contig_fofn = self.get_filepath( "references", "Locus_Contig_Files.fofn" )
        if valid_file( locus_contig_fofn ):
            log.info("Using existing locus contig files and FOFN\n")
            return locus_contig_fofn
        log.info("No locus subread files found, creating...")
        file_prefix = self.get_filepath( "references", "Locus" )
        locus_files = separate_sequences( hla_contigs, contig_loci, file_prefix )
        write_sequence_fofn( locus_files, locus_contig_fofn )
        log.info('Finished separating subreads by locus\n')
        return locus_contig_fofn

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
            print subread_file
            print source
            if not self.to_be_phased( source ):
                continue

            # For sources with a configuration, run AA if output DNE
            output_folder = self.get_filepath( 'phasing', source )
            output_file = os.path.join( output_folder, 'amplicon_analysis.fastq' )
            if os.path.exists( output_file ):
                log.info('AmpliconAnalyzer output detected for "%s", skipping...' % source)
            else:
                log.info('Phasing subreads from "%s"' % source)
                analyzer_args = {'whiteList': subread_file,
                                 'sampleName': '_' + source,
                                 'minLength': self.config(source, 'min_read_length'),
                                 'minReadScore': self.config(source, 'min_read_score'),
                                 'minSnr': self.config(source, 'min_snr'),
                                 #'maxPhasingReads': self.config(source, 'max_phasing_reads'),
                                 'noClustering': self.config(source, 'disable_clustering')}
                analyzer.run( args.input_file, output_folder, analyzer_args )
        log.info('Finished phasing subreads with Amplicon Analysis\n')

    def combine_phasing_results( self ):
        log.info("Looking for the combined output from Amplicon Assembly")
        combined_file = self.get_filepath( 'references', 'amp_analysis_consensus.fastq')
        if valid_file( combined_file ):
            log.info("Using existing combined consensus file\n")
            return combined_file
        log.info("No combined consensus file found, creating...")
        combine_amp_analysis( self.subfolders['phasing'], combined_file )
        log.info('Finished combining AmpliconAnalysis results\n')
        return combined_file

    def run_resequencing( self, baxh5_fofn, renamed_fofn, reference_fofn ):
        log.info("Looking for phased AmpliconAnalysis results")
        print self.smrt_path
        resequencer = Resequencer( self.smrt_path, args.nproc )
        subread_handle = FofnReader( renamed_fofn )
        reference_handle = FofnReader( reference_fofn )
        for subreads, reference in zip(subread_handle, reference_handle):

            # Check if the source of the current file has a configuration
            source = get_file_source( subreads )
            if not self.to_be_resequenced( source ):
                continue

            # For sources with a configuration, run AA if output DNE
            output_folder = self.get_filepath( 'resequencing', source )
            if os.path.exists( output_folder ):
                log.info('Resequencing output detected for "%s", skipping...' % source)
            else:
                log.info('Resequencing HBAR contigs from "%s"' % source)
                resequencer( baxh5_fofn, subreads, reference, output=output_folder )
        log.info('Finished phasing subreads with Amplicon Analysis\n')

    def combine_resequencing_results( self ):
        log.info("Looking for the combined output from Amplicon Assembly")
        combined_file = self.get_filepath( 'references', 'resequencing_consensus.fastq')
        if valid_file( combined_file ):
            log.info("Using existing combined consensus file\n")
            return combined_file
        log.info("No combined consensus file found, creating...")
        combine_resequencing( self.subfolders['resequencing'], combined_file )
        log.info('Finished combining AmpliconAnalysis results\n')
        return combined_file

    def remove_chimeras( self, sequence_file ):
        log.info("Looking for Chimera-filtered AmpliconAnalysis results")
        non_chimera_file = '.'.join( sequence_file.split('.')[:-1] ) + '.good.fastq'
        if valid_file( non_chimera_file ):
            log.info("Using existing Chimera-filtered file\n")
            return non_chimera_file
        log.info("No Chimera-filtered file found, creating...")
        cd = ChimeraDetector()
        cd.run( sequence_file )
        check_output_file( non_chimera_file )
        log.info("Finished removing chimeric sequences\n")
        return non_chimera_file

    def rename_resequencing_results( self, input_file ):
        log.info("Looking for renamed Resequencing consensus file")
        renamed_file = self.get_filepath( 'references', 'resequencing_consensus.renamed.fastq')
        if valid_file( renamed_file ):
            log.info("Using existing renamed consensus file\n")
            return renamed_file
        if not valid_file( input_file ):
            log.info("No valid resequencing output detected, skipping...")
            return input_file
        log.info("No renamed consensus file found, creating...")
        rename_resequencing( input_file, renamed_file )
        check_output_file( renamed_file )
        log.info('Finished combining AmpliconAnalysis results\n')
        return renamed_file

    def combine_phasing_and_reseq( self, good_results, reseq_results ):
        log.info("Looking for combined Phasing and Resequencing results")
        combined_results = self.get_filepath("references", "combined_consensus.fastq")
        if valid_file( combined_results ):
            log.info("Using existing combined consensus file\n")
            return combined_results
        if valid_file( reseq_results ):
            log.info("No combined consensus file found, creating...")
            combine_fastq( [good_results, reseq_results], combined_results )
            check_output_file( combined_results )
        else:
            log.info("No resequencing output to combine, using only Phasing results...")
            copy_file( good_results, combined_results )
        log.info("Finished combining Phasing and Resequencing results\n")
        return combined_results

    def filter_combined_results( self, combined_results ):
        log.info("Looking for filtered consensus results")
        filtered_results = self.get_filepath("references", "filtered_consensus.fastq")
        if valid_file( filtered_results ):
            log.info("Using existing filtered consensus file\n")
            return filtered_results
        log.info("No filtered consensus file found, creating...")
        filter_fastq( combined_results, filtered_results,
                                        min_length=self.min_consensus_length,
                                        min_num_reads=self.min_num_reads)
        check_output_file( filtered_results )
        log.info("Finished filtering combined consensus results\n")
        return filtered_results

    def copy_sequences_for_typing(self, sequence_file):
        log.info("Looking for a sequence file to use in HLA-typing")
        typing_sequences = self.get_filepath('typing', 'consensus_sequences.fastq')
        if valid_file( typing_sequences ):
            log.info("Using existing typing sequences file")
            return typing_sequences
        log.info("No sequence file for typing found, creating...")
        copy_file( sequence_file, typing_sequences )
        log.info("Finished copying sequences for HLA-typing\n")
        return typing_sequences

    def type_hla_sequences(self, sequence_file ):
        log.info('Typing the selected HLA consensus sequences')
        typing = type_sequences( sequence_file,
                                 grouping='locus',
                                 exon_fofn=self.exon_reference,
                                 genomic_reference=self.locus_reference,
                                 cDNA_reference=self.cDNA_reference,
                                 loci=['A', 'B', 'C', 'DQB1', 'DRB1'])
        check_output_file( typing )
        log.info('Finished typing the selected HLA sequences\n')
        return typing

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
                     'alignments', 'phasing', 'typing', 'resequencing']:
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
