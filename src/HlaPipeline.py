#! /usr/bin/env python
import re, os, sys
import shutil
import random
import math
import logging

from subprocess import check_output
from collections import namedtuple, Counter

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbcore.io.GffIO import GffReader, Gff3Record, GffWriter
from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbcore.io.FastqIO import FastqReader, FastqWriter

from pbhla.hbar.HbarTools import HbarRunner
from pbhla.io.BasH5Extractor import BasH5Extractor
from pbhla.io.SamIO import SamReader
from pbhla.io.BlasrIO import parse_blasr
from pbhla.io.GffIO import create_annotation, create_var_annotation
from pbhla.reference.ReferenceDict import ReferenceDict
from pbhla.reference.SequenceSeparator import SequenceSeparator
from pbhla.reference.ReferenceSelector import ReferenceSelector
from pbhla.reference.LocusReference import LocusReference
from pbhla.stats.AmpliconFinder import AmpliconFinder
from pbhla.stats.SubreadStats import SubreadStats
from pbhla.align.MultiSequenceAligner import MSA_aligner
from pbhla.phasing.clusense import Clusense
from pbhla.phasing.SummarizeClusense import combine_clusense_output
from pbhla.resequencing.Resequencer import resequence
from pbhla.external.BlasrRunner import BlasrRunner
from pbhla.external.CdHit import cd_hit_est
from pbhla.external.SmrtAnalysisTools import SmrtAnalysisRunner
from pbhla.fasta.utils import write_fasta, fasta_size, extract_sequence, trim_fasta, combine_fasta
from pbhla.utils import *

__version__ = "0.9.0"
 
file_info = namedtuple('file_info', 'fasta_fn, ref_name, locus')
info = namedtuple('info', 'canonical_pos, feature, codon' )

SMRT_ANALYSIS = "/mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh"
DILUTION = 1.0
MIN_SCORE = 0.8
MIN_LENGTH = 500
NUM_PROC = 4

class HlaPipeline( object ):

    def __init__( self ):
        self.parse_args()
        self.initialize_project()
        self.initialize_logging()
        self.validate_settings()

    def parse_args( self ):
        import argparse
        desc = "A pipeline for performing HLA haplotype sequencing."
        parser = argparse.ArgumentParser( description=desc )

        add = parser.add_argument
        add("input_file", metavar="INPUT",
            help="A BasH5 or FOFN of BasH5s to haplotype.")
        add("reference_file", metavar="REFERENCE", 
            help="FOFN of reference fasta files and their associated loci")
        add("genome", metavar="GENOME", help="Fasta file of the Human Genome")
        add("--output", metavar="DIR", help="Destination folder for process results")
        add("--nproc", type=int, metavar='INT', default=NUM_PROC,
            help="Number of processors to use for parallelization ({0})".format(NUM_PROC))
        add("--min_read_score", type=float, metavar='FLOAT', default=MIN_SCORE,
            help="Only use reads with ReadScore greater than this ({0})".format(MIN_SCORE))
        add("--min_read_length", type=int, metavar='INT', default=MIN_LENGTH,
            help="Only use subreads longer than this ({0})".format(MIN_LENGTH))
        add("--dilution", type=float, metavar='FLOAT', default=DILUTION,
            help="Fraction of subreads to use ({0})".format(DILUTION))
        #add("--MSA", help="FOFN of prealigned MSAs and their associated loci.")
        #add("--region_table", help="Region Table of White-Listed reads to use")
        #add("--phasr-args", nargs='*', default=[''], help="pass these args to phasr.")
        #add("--avoid_phasr", action='store_true',
        #    help="Avoid phasr if reference mapping has identified two likely phases.")
        #add("--annotate", action='store_true', help="Avoid annotation.")
        self.__dict__.update( vars(parser.parse_args()) )

    def initialize_project( self ):
        # If not specified, set the output folder to the input filename
        if self.output is None:
            self.output = self.input_file.split('.')[0]
        create_directory( self.output )
        # Create the various sub-directories
        for d in ['log_files', 'HBAR', 'references', 'subreads', 
                  'alignments', 'phasing', 'phasing_results', 
                  'resequencing', 'stats']:
            sub_dir = os.path.join( self.output, d )
            create_directory( sub_dir )
            setattr(self, d, sub_dir)

    def initialize_logging( self ):
        time_format = "%I:%M:%S"
        log_format = "%(asctime)s %(levelname)s %(filename)s " + \
                     "%(funcName)s %(lineno)d %(message)s"
        out_format = "%(asctime)s %(filename)s %(funcName)s %(message)s"
        self.log = logging.getLogger()
        self.log.setLevel( logging.INFO )
        # Set-up one logger for STDOUT
        h1 = logging.StreamHandler( stream=sys.stdout )
        h1.setLevel( logging.INFO )
        f1 = logging.Formatter( fmt=out_format, datefmt=time_format )
        h1.setFormatter( f1 )
        self.log.addHandler( h1 )
        # Setup a second logger to log to file
        log_file = os.path.join( self.log_files, "HLA_Pipeline.log" )
        h2 = logging.FileHandler( log_file )
        f2 = logging.Formatter( fmt=log_format, datefmt=time_format )
        h2.setFormatter( f2 )
        h2.setLevel( logging.INFO )
        self.log.addHandler( h2 )
    
    def validate_settings( self ):
        # Report the settings with which the pipeline was invoked
        self.log.info( "HLA Pipeline invoked:\n\t{0}\n".format(" ".join(sys.argv)))
        self.input_file = os.path.abspath( self.input_file )
        # Check dilution factors
        if self.dilution <= 0.0 or self.dilution > 1.0:
            msg = "Dilute factor must be between 0 and 1"
            self.log.info( msg )
            raise ValueError( msg )
        # parse phasr args
        #self.phasr_argstring = ''
        #for argument in self.phasr_args:
        #    if ':' in argument:
        #        param, value = argument.split(":")
        #        self.phasr_argstring += '--%s %s ' % (param, value)
        #    elif 'output' in argument or 'cname' in argument:
        #        pass    
        #    else:
        #        self.phasr_argstring += '--%s ' % argument
        # Initialize the SmrtAnalysisTools
        self.smrt_analysis = SmrtAnalysisRunner( SMRT_ANALYSIS, self.log_files, self.nproc )

    def getVersion(self):
        return __version__

    def __call__( self ):
        # First we assemble the supplied data via HGAP / HBAR
        self.assemble_contigs()
        # Second we extract the subreads ourselves and map them onto
        #     the created contigs
        self.extract_subreads()
        self.align_subreads_to_contigs()
        self.create_subread_contig_dict()
        # Third, we align the contigs to the Human Genome ...
        self.align_contigs_to_genome()
        self.create_contig_genome_dict()
        # ... and the locus reference sequences ...
        self.create_ref_locus_dict()
        self.create_locus_reference()
        self.align_contigs_to_reference()
        self.create_contig_reference_dict()
        self.create_contig_locus_dict()
        # ... in order to separate the on-target from off-target hits
        self.find_on_target_contigs()
        self.separate_off_target_contigs()
        # Fourth, we combine the Subread and Contig dicts and use them
        #     to separate out off-target subreads
        self.find_on_target_subreads()
        self.separate_off_target_subreads()
        # Next we create re-align just the HLA data for summarizing
        self.remove_redundant_contigs()
        self.realign_hla_subreads()
        # Next we summarize our pre-phasing coverage of the contigs
        self.summarize_aligned_hla_subreads()
        self.separate_contigs()
        self.create_subread_selected_contig_dict()
        self.separate_subreads_by_contig()
        self.phase_reads_with_clusense()
        self.summarize_clusense()
        self.output_phased_contigs()
        self.resequence_contigs()
        self.separate_resequenced_contigs()
        self.realign_subreads_to_resequenced()
        self.align_resequenced_to_reference()
        self.create_resequenced_dict()
        self.create_resequenced_locus_dict()
        self.create_subread_resequenced_contig_dict()
        self.separate_subreads_by_resequenced_contig()

    def check_output_file(self, filepath):
        try:
            assert os.path.isfile( filepath )
        except:
            msg = 'Expected output file not found! "{0}"'.format(filepath)
            self.log.info( msg )
            raise IOError( msg )

    def assemble_contigs( self ):
        self.log.info('Beginning the extraction of subiread data')
        contig_output = os.path.join( self.HBAR,
                                      "3-CA",
                                      "9-terminator",
                                      "asm.utg.fasta" )
        self.contig_file = os.path.join( self.references,
                                         "all_contigs.fasta")
        if os.path.isfile( self.contig_file ):
            self.log.info('Found existing contig file "{0}"'.format(self.contig_file))
            self.log.info('Skipping HBAR assembly step\n')
            return
        # Run HGAP
        self.log.info('No contig file found, initializing new HbarRunner')
        hbar = HbarRunner( self.input_file, self.HBAR )
        hbar()
        # Copy the contig file to a more convenient location
        shutil.copy( contig_output, self.contig_file )
        self.check_output_file( self.contig_file )
        self.log.info('Finished the assembly of subread data\n')

    def extract_subreads( self ):
        self.log.info('Beginning the extraction of subread data')
        # Dump all valid reads from the above files
        self.subread_file = os.path.join( self.subreads, "all_subreads.fasta" )
        if os.path.isfile( self.subread_file ):
            self.log.info('Found existing subread file "{0}"'.format(self.subread_file))
            self.log.info('Skipping subread extraction step\n')
            return
        self.log.info('Extracting subreads from input files...')
        BasH5Extractor( self.input_file, 
                        self.subread_file,
                        min_length=self.min_read_length,
                        min_score=self.min_read_score,
                        dilution=self.dilution)
        self.check_output_file( self.subread_file )
        self.log.info('Finished the extraction of subread data\n')

    def align_subreads_to_contigs(self):
        self.log.info('"Aligning all subreads to the HBAR contigs')
        self.subreads_to_contigs = os.path.join(self.alignments, 'subreads_to_contigs.m1')   
        if os.path.isfile( self.subreads_to_contigs ):
            self.log.info('Found existing alignment file "{0}"'.format(self.subreads_to_contigs))
            self.log.info("Skipping alignment step...\n")
            return
        subread_count = fasta_size( self.subread_file )
        contig_count = fasta_size( self.contig_file )
        # Run BLASR
        self.log.info("Aligning {0} contigs to the {1} contigs".format(subread_count, contig_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': 30}
        BlasrRunner( self.subread_file, self.contig_file, self.subreads_to_contigs, blasr_args )
        # Check and save the output
        self.check_output_file( self.subreads_to_contigs )
        self.log.info("Finished aligning subreads to the created contigs\n")

    def create_subread_contig_dict(self):
        self.log.info("Converting the Subreads-to-Contigs alignment to a Dictionary")
        self.subread_contig_dict = ReferenceDict( self.subreads_to_contigs )
        self.log.info("Finished convereting the data to a Dictionary\n")

    def align_contigs_to_genome(self):
        self.log.info('"Aligning all contigs to the Human Genome')
        self.contigs_to_genome = os.path.join(self.alignments, 'contigs_to_genome.m1')   
        if os.path.isfile( self.contigs_to_genome ):
            self.log.info('Found existing alignment file "{0}"'.format(self.contigs_to_genome))
            self.log.info("Skipping alignment step...\n")
            return
        sa_file = self.genome + '.sa'
        contig_count = fasta_size( self.contig_file )
        # Run BLASR
        self.log.info("Aligning {0} contigs to the genomic reference".format(contig_count))
        blasr_args = {'nproc': self.nproc,
                      'sa': sa_file,
                      'bestn': 1,
                      'nCandidates': 10}
        BlasrRunner( self.contig_file, self.genome, self.contigs_to_genome, blasr_args )
        # Check and save the output
        self.check_output_file( output_path )
        self.log.info("Finished aligning contigs to the genomic reference\n")

    def create_contig_genome_dict(self):
        self.log.info("Converting the Contigs-to-Genome alignment to a Dictionary")
        self.contig_genome_dict = ReferenceDict( self.contigs_to_genome )
        self.log.info("Finished convereting the data to a Dictionary\n")

    def create_ref_locus_dict(self):
        self.log.info('"Creating reference dictonary from "{0}"'.format(self.reference_file))
        self.reference_locus_dict = ReferenceDict( self.reference_file )
        self.log.info("Finished creating locus reference dict...\n")

    def create_locus_reference(self):
        self.log.info("Creating locus reference fasta file")
        self.reference_seqs = os.path.join(self.references, "locus_references.fasta")
        # If no locus key and Choose_Ref is None, read the locus from the regular reference
        if os.path.isfile( self.reference_seqs ):
            self.log.info('Found existing locus reference sequences "{0}"'.format(self.reference_seqs))
            self.log.info('Skipping reference extraction step\n')
            return
        self.log.info("No locus reference sequences found, creating one...")
        LocusReference(self.reference_file, self.reference_seqs)
        self.check_output_file( self.reference_seqs )
        self.log.info("Finished creating locus reference\n")
    
    def align_contigs_to_reference(self):
        self.log.info('"Aligning all HLA contigs to the Reference sequences')
        self.contigs_to_reference = os.path.join(self.alignments, 'contigs_to_reference.m1')   
        if os.path.isfile( self.contigs_to_reference ):
            self.log.info('Found existing alignment file "{0}"'.format(self.contigs_to_reference))
            self.log.info("Skipping alignment step...\n")
            return
        contig_count = fasta_size( self.contig_file )
        reference_count = fasta_size( self.reference_seqs )
        # Run BLASR
        self.log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': reference_count}
        BlasrRunner( self.contig_file, self.reference_seqs, self.contigs_to_reference, blasr_args )
        # Check and save the output
        self.check_output_file( self.contigs_to_reference )
        self.log.info("Finished aligning contigs to the HLA reference set\n")

    def create_contig_reference_dict(self):
        self.log.info("Converting the Contigs-to-Reference alignment to a Dictionary")
        self.contig_reference_dict = ReferenceDict( self.contigs_to_reference )
        self.log.info("Finished convereting the data to a Dictionary\n")

    def create_contig_locus_dict(self):
        self.log.info("Converting the Contigs-to-Reference dict and Reference-to-Locus")
        self.log.info("    dict into a combined Contig-to-Locus Dictionary")
        self.contig_locus_dict = cross_ref_dict( self.contig_reference_dict, self.reference_locus_dict )
        self.log.info("Finished convereting the data to a Dictionary\n")

    def find_on_target_contigs(self):
        self.log.info('Identifying on-target contigs')
        self.on_target_contig_ids = [c for c in self.contig_locus_dict
                                       if self.contig_genome_dict[c] == 'chr6']
        self.log.info('Finished identifying on-target contigs\n')

    def separate_off_target_contigs(self):
        self.log.info('Separating off-target contigs by genomic alignment')
        self.hla_contigs = os.path.join(self.references, 'hla_contigs.fasta')
        self.non_hla_contigs = os.path.join(self.references, 'non_hla_contigs.fasta')
        if os.path.isfile( self.hla_contigs ) and os.path.isfile( self.non_hla_contigs ):
            self.log.info('Found existing sequence file "{0}"'.format(self.hla_contigs))
            self.log.info('Found existing sequence file "{0}"'.format(self.non_hla_contigs))
            self.log.info("Skipping separation step...\n")
            return
        self.log.info('No separated contig files found, initializing separator')
        separator = SequenceSeparator( self.contig_file, selected=self.on_target_contig_ids )
        separator.write('selected', self.hla_contigs)
        separator.write('not_selected', self.non_hla_contigs)
        self.log.info('Finished separating off-target contigs\n')
    
    def find_on_target_subreads(self):
        self.log.info('Identifying on-target subreads')
        self.on_target_subreads = [s for s in self.subread_contig_dict
                                   if self.subread_contig_dict[s] in self.on_target_contig_ids]
        self.log.info('Finished identifying on-target subreads\n')

    def separate_off_target_subreads(self):
        self.log.info('Separating off-target subreads by contig alignment')
        self.hla_subreads = os.path.join(self.subreads, 'hla_subreads.fasta')
        self.non_hla_subreads = os.path.join(self.subreads, 'non_hla_subreads.fasta')
        if os.path.isfile( self.hla_subreads ) and os.path.isfile( self.non_hla_subreads ):
            self.log.info('Found existing sequence file "{0}"'.format(self.hla_subreads))
            self.log.info('Found existing sequence file "{0}"'.format(self.non_hla_subreads))
            self.log.info("Skipping separation step...\n")
            return
        self.log.info('No separated contig files found, initializing separator')
        separator = SequenceSeparator( self.subread_file, selected=self.on_target_subreads )
        separator.write('selected', self.hla_subreads)
        separator.write('not_selected', self.non_hla_subreads)
        self.log.info('Finished separating off-target subreads\n')

    def remove_redundant_contigs(self):
        self.log.info("Removing redundant contigs from the HLA contig file")
        self.selected_contigs = os.path.join(self.references, 'selected_contigs.fasta')
        if os.path.isfile( self.selected_contigs ):
            self.log.info('Found existing non-redundant contig file "{0}"'.format(self.selected_contigs))
            self.log.info("Skipping redundant contig filtering step...\n")
            return
        self.log.info('File of selected contigs not found, initializing Cd-Hit-Est')
        cd_hit_est( self.hla_contigs, self.selected_contigs )
        self.log.info('Finished selecting non-redundant contigs\n')

    def realign_hla_subreads(self):
        self.log.info("Re-aligning HLA subreads to selected references")
        self.hla_alignment = os.path.join( self.alignments, "hla_subreads_to_contigs.sam" )
        if os.path.isfile( self.hla_alignment ):
            self.log.info('Found existing SAM file "{0}"'.format(self.hla_alignment))
            self.log.info("Skipping realignment step...\n")
            return
        query_count = fasta_size( self.hla_subreads )
        ref_count = fasta_size( self.selected_contigs )
        self.log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': ref_count}
        BlasrRunner(self.hla_subreads, self.selected_contigs, self.hla_alignment, blasr_args)
        self.check_output_file( self.hla_alignment )
        self.log.info("Finished realigning HLA subreads to HLA contigs\n")

    def separate_contigs(self):
        self.log.info("Separating remaining contigs into individual files")
        contig_fofn = os.path.join( self.subreads, "contig_files.txt" )
        if os.path.isfile( contig_fofn ):
            self.log.info('Found existing Contig File List "{0}"'.format(contig_fofn))
            self.log.info("Skipping contig separation step...\n")
            self.contig_files = read_fofn( contig_fofn )
            return
        separator = SequenceSeparator( self.selected_contigs )
        contig_prefix = os.path.join(self.references, 'Contig')
        self.contig_files = separator.write_all( contig_prefix )
        write_fofn( self.contig_files, contig_fofn)
        self.log.info('Finished separating contigs')

    def create_subread_selected_contig_dict(self):
        self.log.info("Converting the Subreads-to-Contigs alignment to a Dictionary")
        self.subread_selected_contig_dict = ReferenceDict( self.hla_alignment )
        self.log.info("Finished convereting the data to a Dictionary\n")

    def separate_subreads_by_contig(self):
        self.log.info("Separating subreads by aligned contig")
        subread_fofn = os.path.join( self.subreads, "subread_files.txt" )
        if os.path.isfile( subread_fofn ):
            self.log.info('Found existing Subread File List "{0}"'.format(subread_fofn))
            self.log.info("Skipping subread separation step...\n")
            self.subread_files = read_fofn( subread_fofn )
            return
        separator = SequenceSeparator( self.hla_subreads, 
                                       reference_dict=self.subread_selected_contig_dict )
        subread_prefix = os.path.join(self.subreads, 'Subreads')
        self.subread_files = separator.write_all( subread_prefix )
        self.subread_files = [fn for fn in self.subread_files if not fn.endswith('Unmapped.fasta')]
        write_fofn( self.subread_files, subread_fofn )
        self.log.info('Finished separating subreads by contig')

    def phase_reads_with_clusense(self):
        self.log.info("Phasing subreads with Clusense")
        for subread_fn, contig_fn in zip(self.subread_files, self.contig_files):
            folder_name = os.path.basename(contig_fn).split('.')[0]
            output_folder = os.path.join(self.phasing, folder_name)
            if os.path.exists( output_folder ):
                self.log.info('"Clusense output detected for "{0}", skipping...'.format(folder_name))
            else:
                self.log.info("Phasing subreads for {0}".format(folder_name))
                Clusense(subread_fn, contig_fn, output_folder)
        self.log.info('Finished Phasing subreads with Clusense')

    def summarize_clusense(self):
        self.log.info("Summarizing the output from Clusense")
        output_folder = os.path.join(self.phasing_results, 'Clusense_Results')
        self.clusense_cns_files, self.clusense_read_files = combine_clusense_output(self.phasing, output_folder)
        self.log.info('Finished Phasing subreads with Clusense')

    def output_phased_contigs(self):
        self.log.info("Creating a combined reference of non-HLA and phased-HLA contigs")
        self.phased_contigs = os.path.join(self.references, 'phased_contigs.fasta')
        fasta_files = list(self.clusense_cns_files) + [self.non_hla_contigs]
        combine_fasta(fasta_files, self.phased_contigs)
        self.log.info('Finished creating combined reference')

    def resequence_contigs(self):
        self.log.info("Resequencing phased contigs")
        self.resequenced_contigs = os.path.join(self.references, 'resequenced_contigs.fasta')
        if os.path.exists( self.resequenced_contigs ):
            self.log.info('Existing Quiver consensus output found, skipping...')
            return
        consensus_file = resequence( self.smrt_analysis, self.input_file, self.phased_contigs, self.resequencing)
        shutil.copy( consensus_file, self.resequenced_contigs )
        self.log.info('Finished the contig resequencing step')

    def separate_resequenced_contigs(self):
        self.log.info("Resequencing phased contigs")
        self.resequenced_hla = os.path.join(self.references, 'resequenced_hla.fasta')
        self.resequenced_non_hla = os.path.join(self.references, 'resequenced_non_hla.fasta')
        if os.path.isfile( self.resequenced_hla ):
            self.log.info('Found existing file of resequenced HLA contigs')
            self.log.info('Skipping resequenced contig separation step...\n')
            return
        with FastaWriter(self.resequenced_hla) as hla:
            with FastaWriter(self.resequenced_non_hla) as non_hla:
                for record in FastaReader( self.resequenced_contigs ):
                    if record.name.startswith('Contig_'):
                        hla.writeRecord( record )
                    else:
                        non_hla.writeRecord( record )
        self.log.info('Finished separating the resequenced contigs')

    def align_resequenced_to_reference(self):
        self.log.info('"Aligning all resequenced HLA contigs to the Reference sequences')
        self.resequenced_to_reference = os.path.join(self.alignments, 'resequenced_to_reference.m1')   
        if os.path.isfile( self.resequenced_to_reference ):
            self.log.info('Found existing alignment file "{0}"'.format(self.resequenced_to_reference))
            self.log.info("Skipping alignment step...\n")
            return
        contig_count = fasta_size( self.resequenced_hla )
        reference_count = fasta_size( self.reference_seqs )
        # Run BLASR
        self.log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': reference_count}
        BlasrRunner( self.resequenced_hla, self.reference_seqs, self.resequenced_to_reference, blasr_args )
        # Check and save the output
        self.check_output_file( self.resequenced_to_reference )
        self.log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def create_resequenced_dict(self):
        self.log.info("Converting the Resequenced-to-Reference alignment to a Dictionary")
        self.resequenced_reference_dict = ReferenceDict( self.resequenced_to_reference )
        self.log.info("Finished convereting the data to a Dictionary\n")

    def create_resequenced_locus_dict(self):
        self.log.info("Converting the Contigs-to-Reference dict and Reference-to-Locus")
        self.log.info("    dict into a combined Contig-to-Locus Dictionary")
        self.resequenced_locus_dict = cross_ref_dict( self.resequenced_reference_dict, self.reference_locus_dict )
        self.log.info("Finished convereting the data to a Dictionary\n")

    def realign_subreads_to_resequenced(self):
        self.log.info("Re-aligning HLA subreads to resequenced references")
        self.resequenced_alignment = os.path.join( self.alignments, "hla_subreads_to_resequenced.sam" )
        if os.path.isfile( self.resequenced_alignment ):
            self.log.info('Found existing SAM file "{0}"'.format(self.resequenced_alignment))
            self.log.info("Skipping realignment step...\n")
            return
        query_count = fasta_size( self.hla_subreads )
        ref_count = fasta_size( self.resequenced_hla )
        self.log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': ref_count}
        BlasrRunner(self.hla_subreads, self.resequenced_hla, self.resequenced_alignment, blasr_args)
        self.check_output_file( self.resequenced_alignment )
        self.log.info("Finished realigning HLA subreads to resequenced contigs\n")

    def create_subread_resequenced_contig_dict(self):
        self.log.info("Converting the Subreads-to-Resequenced-Contigs alignment to a Dictionary")
        self.subread_resequenced_contig_dict = ReferenceDict( self.resequenced_alignment )
        self.log.info("Finished convereting the data to a Dictionary\n")

    def separate_subreads_by_resequenced_contig(self):
        self.log.info("Separating subreads by resequenced contig")
        self.subread_fofn = os.path.join( self.subreads, "resequenced_subread_files.txt" )
        if os.path.isfile( subread_fofn ):
            self.log.info('Found existing Subread File List "{0}"'.format(self.subread_fofn))
            self.log.info("Skipping subread separation step...\n")
            self.subread_files = read_fofn( self.subread_fofn )
            return
        separator = SequenceSeparator( self.hla_subreads, 
                                       reference_dict=self.subread_resequenced_contig_dict )
        subread_prefix = os.path.join(self.subreads, 'Resequenced')
        self.subread_files = separator.write_all( subread_prefix )
        self.subread_files = [fn for fn in self.subread_files if not fn.endswith('Unmapped.fasta')]
        write_fofn( self.subread_files, self.subread_fofn )
        self.log.info('Finished separating subreads by contig')

    def pick_final_contigs(self):
        self.log.info("Picking contigs from the resequencing and realignment results")
        ContigPicker( self.resequenced_hla, self.subread_fofn, self.resequenced_locus_dict )
        self.log.info("Finished picking contigs")

    def summarize_aligned_hla_subreads(self):
        self.log.info("Summarizing coverage of the HLA contigs by the HLA subreads")
        self.locus_stats = os.path.join(self.stats, "locus_statistics.csv")
        self.reference_stats = os.path.join(self.stats, "reference_statistics.csv")
        ### will go through sam file, sort out the raw reads, and tabulate statistics while we do it
        if os.path.isfile( self.locus_stats ) and \
           os.path.isfile( self.reference_stats ):
            self.log.info('Found existing Locus Statistics at "{0}"'.format(self.locus_stats))
            self.log.info('Found existing Reference Statistics at "{0}"'.format(self.reference_stats))
            self.log.info("Skipping alignment summary step...\n")
            return
        stats = SubreadStats( self.hla_contigs, self.contig_locus_dict )
        stats.add_sam_file( self.hla_alignment )
        ### finally write out the subread statistics
        stats.write( self.locus_stats, 'locus' )
        stats.write( self.reference_stats, 'reference' )
        self.log.info("Finished summarizing coverage of the HLA contigs")

    def find_amplicons(self):
        self.log.info("Finding amplicon locations within each reference")
        self.amplicon_summary = self.stats_dir + "/amplicon_summary.csv"
        if os.path.isfile( self.amplicon_summary ):
            self.log.info('Found existing Amplicon Summary "{0}" )'.format(self.amplicon_summary))
            self.log.info("Skipping amplicon finding step...\n")
            return
        AmpliconFinder(self.sam_file, self.amplicon_summary, self.locus_dict)
        self.log.info("Finished finding amplicon locations\n")

    def phase_subreads( self ):
        self.log.info("Starting the sub-read phasing process")
        # First 
        subread_files = []
        with open( self.subread_files ) as handle:
            for line in handle:
                info_tuple = file_info._make( line.strip().split() )
                subread_files.append( info_tuple )
        ### run through the subreads files to figure out which loci should actually be run through phasr
        ### if there are already two subread files for a single locus, theres little point in running through phasr
        ### as the two alleles are already well separated
        loci_count = Counter([ x.locus for x in subread_files ])
        to_be_phased={}
        for item in loci_count.iteritems():
            if item[1] == 1:
                to_be_phased[item[0]] = False
            else:
                to_be_phased[item[0]] = True
        
        seq_to_locus={}
        backbone_fasta = self.args.proj+"/phased/backbone.fa"
        self.hap_con_fasta = self.args.proj + "/phased/haplotype_consensus_sequences.fasta"
        for item in subread_files:
            ref_name = item.ref_name
            locus = item.locus
            fasta_file = item.fasta_fn
            ### if these are the unmapped reads we are going to run them through HGAp
            if locus == "unmapped":
                if self.args.HGAP:
                    pass
                    #phase_unmapped_reads()
                continue
            
            phasr_output = self.args.proj+"/phased/"+ref_name+"_haplotype_consensus.fa"

            ### extract the reference sequence that was used to "cluster" this group of reads
            sequence = extract_sequence(self.reference_sequences, [ ref_name ] )
            write_fasta(sequence, backbone_fasta, "w")
            
            try:
                assert os.path.isfile( fasta_file )
            except:
                msg = "Could not open ( %s ). Skipping." % fasta_file
                self.log.info( msg )
                raise IOError( msg )

            if os.path.isfile(phasr_output):
                for record in FastaReader(phasr_output):
                    seq_to_locus[record.name] = locus
                self.log.info("Found phasr output for ( %s ) in ( %s ). Skipping." % (fasta_file, phasr_output))
                continue
            
            ### if this locus is not to be phased, we can just set max recursion level to 0, and phasr will build regular consensus
            if self.args.avoid_phasr and not to_be_phased[locus]:
                self.log.info('Running Gcon instead of Phasr for "%s"' % fasta_file)
                argstring = re.sub('--max_recursion_level [0-9]',
                                   '--max_recursion_level 0',
                                   self.phasr_argstring)
                argstring += ' --max_recursion_level 0'
            else:
                self.log.info('Running  Phasr on "%s"' % fasta_file)
                argstring = self.phasr_argstring
            phasr_cmd = "phasr %s --ref %s\
                --output %s \
                --log \
                %s" % ( fasta_file, backbone_fasta, phasr_output, self.phasr_argstring)
            ### run phasr
            self.log.info("phasr invoked: %s" % phasr_cmd)
            check_output(". /home/UNIXHOME/jquinn/HGAP_env/bin/activate; %s" % (phasr_cmd),  
                         executable='/bin/bash', shell=True)
            if os.path.isfile( phasr_output ):
                output_count = fasta_size(phasr_output)
                self.log.info('Phasr output %s seqs to "%s"' % ( output_count, phasr_output) )
                ### append results to multifasta
                runbash( "cat %s >> %s" % (phasr_output, self.hap_con_fasta) )
                for record in FastaReader(phasr_output):
                    seq_to_locus[record.name] = locus
                ### create record of files created
                with open(self.args.proj+"/phased/phasr_output_files.txt", "a") as of:
                    print >>of, "%s %s %s" % (phasr_output, ref_name, locus)
        # Finally, we write a summary of our Phasr outputs to file
        with open(self.args.proj+"/phased/phasr_output_seqs.txt", "w") as of:
            for item in seq_to_locus.iteritems():
                print >>of, "%s %s" % (item[0], item[1])
        phasr_output_count = fasta_size( self.hap_con_fasta )
        self.log.info("Phasr created %s sequence(s) total" % phasr_output_count )
        self.log.info("Phasing complete.\n")

    def annotate(self):     
        if not self.args.annotate:
            return
        try:
            os.mkdir(self.args.proj+"/annotate")
        except:
            pass
        self.log.info("Annotation begun")
        self.locus_dict={}       
        nseqs=0
        with open(self.args.proj+"/phased/phasr_output_seqs.txt", "r") as f:
            for line in f:
                self.locus_dict[line.strip().split()[0]] = line.strip().split()[1]
                nseqs+=1
        self.log.info("( %s ) sequences to annotate." % (nseqs) )
        MSA_fn_dict={}; 
        MSA_cDNA_fn_dict={}
        MSA_info_fn_dict={}; 
        MSA_cDNA_info_fn_dict={}
        tmp_fn = self.args.proj+"/annotate/tmp.fasta"
        with open(self.args.MSA, "r") as f:
            for line in f:
                line = line.split()
                MSA_fn_dict[line[2]] = line[0]+".afa"
                MSA_info_fn_dict[line[2]] = line[0]+".info"
                MSA_cDNA_info_fn_dict[line[2]] = line[1]+".info"
                MSA_cDNA_fn_dict[line[2]] = line[1]+".afa"
        gDNA_var_writer = GffWriter(self.args.proj+"/annotate/best_gDNA_ref_match_variants.gff")
        gDNA_var_writer.writeMetaData('pacbio-variant-version', '1.4')      
        cDNA_var_writer = GffWriter(self.args.proj+"/annotate/best_cDNA_ref_match_variants.gff")
        cDNA_var_writer.writeMetaData('pacbio-variant-version', '1.4')
        gDNA_annot_writer = GffWriter(self.args.proj+"/annotate/gDNA.gff")
        cDNA_annot_writer = GffWriter(self.args.proj+"/annotate/cDNA.gff")
        feature_writer = GffWriter(self.args.proj+"/annotate/features.gff")
        gDNA_calls={}
        cDNA_calls={}
        if os.path.isfile(self.args.proj+"/annotate/gDNA_allele_calls.txt"):
            os.remove(self.args.proj+"/annotate/gDNA_allele_calls.txt")
        if os.path.isfile(self.args.proj+"/annotate/resequenced_hap_con_cDNA.fasta"):
            os.remove(self.args.proj+"/annotate/resequenced_hap_con_cDNA.fasta")
        if os.path.isfile(self.args.proj+"/annotate/cDNA_allele_calls.txt"):
            os.remove(self.args.proj+"/annotate/cDNA_allele_calls.txt")
        f = FastqReader(self.args.proj+"/reseq/resequenced_hap_con.fastq")
        for r in f:
            write_fasta([r], tmp_fn, "w") 
            name = r.name.split("|")[0]
            locus = self.locus_dict[name]
            output = self.args.proj+"/annotate/"+name+".afa"
            self.log.info("Processing ( %s ) from locus ( %s )." % (name, locus) )

            ### read in profile features
            MSA_info_fn = MSA_info_fn_dict[locus]
            MSA_info = {}
            with open(MSA_info_fn, "r") as f2:
                for line in f2:
                    line = line.strip().split()
                    MSA_info[int(line[0])] = info._make( ( line[1], line[2], 0 ) )

            ### read profile fn
            MSA_fn = MSA_fn_dict[locus]

            ### align sequence to profile for this locus
            if not os.path.isfile(output):
                try:
                    muscle_output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
                        muscle -profile -in1 %s -in2 %s -out %s 2> /dev/null" % (MSA_fn, tmp_fn, output),
                        executable='/bin/bash', shell=True)
                except:
                    self.log.info("MSA failed for ( %s )" % (name) )
                    os.remove(tmp_fn)
                    continue
            ### read out a reference sequence from the old profile
            f2 = FastaReader(MSA_fn)
            letters_max = 0
            for r2 in f2:
                num_letters = ( r2.sequence.count('A') + r2.sequence.count('G') + r2.sequence.count('C') + r2.sequence.count('T') )
                if num_letters > letters_max:
                    letters_max = num_letters
            f2 = FastaReader(MSA_fn)
            for r2 in f2:
                num_letters = ( r2.sequence.count('A') + r2.sequence.count('G') + r2.sequence.count('C') + r2.sequence.count('T') )
                if num_letters == letters_max:
                    comparison_seq = r2.sequence
                    comparison_seq_name = r2.name
                    break

            ### read out the same reference sequence from the NEW profile
            ### and read out our consensus sequence
            f2 = FastaReader(output)
            for r2 in f2:
                if r2.name == comparison_seq_name:
                    new_comparison_seq = r2.sequence
                    new_comparison_seq_name = r2.name
                elif r2.name == r.name: 
                    consensus_seq = r2.sequence

            ### read coverage info from the cov gff for this sequence
            f2 = GffReader(self.args.proj+"/reseq/coverage_1.gff")
            coverage_map={}
            for record in f2:
                if record.seqid == name:
                    for i in xrange(int(record.start), int(record.end)+1):
                        coverage_map[i] = record.getAttrVal('cov2') 

            ### read Qv info from the fastq for this sequence
            quality_map = {}
            i = 1
            for qual in r.quality:
                quality_map[i] = qual
                i += 1

            ### create annotation gff3 for this sequence
            ### comparing the same ref between the new and old profile allows us to keep track of the "canonical" coordinates
            gff3_annotation = create_annotation(comparison_seq, new_comparison_seq, consensus_seq, coverage_map, MSA_info, r.name, quality_map)     
            last_feature = None; feature_regions = {}
            for record in gff3_annotation:
                gDNA_annot_writer.writeRecord(record)       
                try:
                    feature_regions[record.getAttrVal('feature')].append(record.getAttrVal('sequence_position'))
                except KeyError:
                    feature_regions[record.getAttrVal('feature')] = [int(record.getAttrVal('sequence_position'))]
            for item in feature_regions.iteritems():
                feature = item[0]
                if 'exon' not in feature: continue
                bases = sorted(item[1])     
                record = Gff3Record(r.name) 
                record.type='region'
                record.start=bases[0]
                record.end=bases[-1]
                record.put('feature', feature )
                feature_writer.writeRecord(record)
                
            ### find best match among all profiles for our consensus sequence       
            ### write it out to a file
            best_gDNA_match, best_gDNA_match_score, gDNA_var_map = MSA_aligner(consensus_seq, output, r.name)           
            with open(self.args.proj+"/annotate/gDNA_allele_calls.txt", "a") as of:
                print>>of, "%s %s %s" % (r.name, best_gDNA_match, best_gDNA_match_score)    
            gDNA_calls[r.name] = locus
            
            ### create variant gff3 to document differences between consensus sequence and best matching reference
            gff3_var_annotation = create_var_annotation(consensus_seq, gDNA_var_map, r.name)    
            for record in gff3_var_annotation:
                gDNA_var_writer.writeRecord(record)

            ### create artificial cDNA sequence from the gDNA sequence
            ### create a dict between gDNA and cDNA position so that we can read in quality and coverage info to the cDNA gff
            ### and so we can keep track of where our artificial cDNA sequence fits in to the real sequence
            cDNA_consensus_sequence=''
            cDNA_base_counter=1
            gDNA_pos_to_cDNA_pos={}
            for i in xrange(len(gff3_annotation)):
                if gff3_annotation[i].type == 'base':       
                    if gff3_annotation[i].getAttrVal('isexon'):
                        cDNA_consensus_sequence+=gff3_annotation[i].getAttrVal('basecall')
                        gDNA_pos_to_cDNA_pos[gff3_annotation[i].getAttrVal('sequence_position')] = cDNA_base_counter
                        cDNA_base_counter+=1
            cDNA_coverage_map={}
            cDNA_pos_to_gDNA_pos={}
            for item in gDNA_pos_to_cDNA_pos.iteritems():
                cDNA_pos_to_gDNA_pos[int(item[1])] = item[0]
            
            for item in coverage_map.iteritems():
                try:
                    cDNA_coverage_map[gDNA_pos_to_cDNA_pos[item[0]]] = item[1]
                except:
                    pass

            cDNA_quality_map={}
            for item in quality_map.iteritems():
                try:
                    cDNA_quality_map[gDNA_pos_to_cDNA_pos[item[0]]] = item[1]
                except:
                    pass

            with open(self.args.proj+"/annotate/resequenced_hap_con_cDNA.fasta", "a") as of:
                print>>of, ">"+r.name
                print>>of, cDNA_consensus_sequence

            ### now that we have established the cDNA sequence ... we start all over again with a cDNA MSA!
            ### get cDNA profile fn
            MSA_cDNA_info_fn = MSA_cDNA_info_fn_dict[locus]

            ### read in cDNA profile features
            MSA_cDNA_info = {}
            with open(MSA_cDNA_info_fn, "r") as f2:
                for line in f2:
                    line = line.strip().split()
                    MSA_cDNA_info[int(line[0])] = info._make( ( line[1], line[2], line[3] ) )
            MSA_cDNA_fn = MSA_cDNA_fn_dict[locus]
            with open(tmp_fn, "w") as of:
                print>>of, ">"+r.name       
                print>>of, cDNA_consensus_sequence

            ### add artificial cDNA sequence to the cDNA profile
            output = self.args.proj+"/annotate/"+name+"_cDNA.afa"
            if not os.path.isfile(output):
                try:
                    muscle_output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
                        muscle -profile -in1 %s -in2 %s -out %s 2> /dev/null" % (MSA_cDNA_fn, tmp_fn, output),
                        executable='/bin/bash', shell=True)
                except:
                    self.log.info("MSA failed for ( %s )" % (name+"_cDNA") )
                    os.remove(tmp_fn)
                    continue

            ### read out a reference sequence from the profile
            f2 = FastaReader(MSA_cDNA_fn)
            letters_max = 0
            for r2 in f2:
                num_letters = ( r2.sequence.count('A') + r2.sequence.count('G') + r2.sequence.count('C') + r2.sequence.count('T') )
                if num_letters > letters_max:
                    letters_max = num_letters
            f2 = FastaReader(MSA_cDNA_fn)
            for r2 in f2:
                num_letters = ( r2.sequence.count('A') + r2.sequence.count('G') + r2.sequence.count('C') + r2.sequence.count('T') )
                if num_letters == letters_max:
                    comparison_seq = r2.sequence
                    comparison_seq_name = r2.name
                    break

            ### read out the same reference sequence from the NEW profile
            ### as well as our artificial cDNA consensus sequence
            f2 = FastaReader(output)
            for r2 in f2:
                if r2.name == comparison_seq_name:
                    new_comparison_seq = r2.sequence
                    new_comparison_seq_name = r2.name
                elif r2.name == r.name: 
                    consensus_seq = r2.sequence

            ### create annotation gff3 for this sequence
            ### comparing the same ref between the new and old profile allows us to keep track of the "canonical" coordinates
            gff3_annotation = create_annotation(comparison_seq, new_comparison_seq, consensus_seq, cDNA_coverage_map, MSA_cDNA_info, r.name, cDNA_quality_map, coord_dict = cDNA_pos_to_gDNA_pos )  
            for record in gff3_annotation:
                cDNA_annot_writer.writeRecord(record)
            ## find best match among all profiles for our consensus sequence    
            best_cDNA_match, best_cDNA_match_score, cDNA_var_map = MSA_aligner(consensus_seq, output, r.name)           
            with open(self.args.proj+"/annotate/cDNA_allele_calls.txt", "a") as of:
                print>>of, "%s %s %s" % (r.name, best_cDNA_match, best_cDNA_match_score)
            cDNA_calls[r.name] = locus

            ### create variant gff3 to document differences between consensus sequence and best matching reference
            gff3_var_annotation = create_var_annotation(consensus_seq, cDNA_var_map, r.name, coord_dict = cDNA_pos_to_gDNA_pos )    
            for record in gff3_var_annotation:
                cDNA_var_writer.writeRecord(record)

            os.remove(tmp_fn)
        ### now we create a special gff that will allow us to do phase highlighting
        ### TODO: need to have a locus dict
        locus_counts = Counter(self.locus_dict.values())
        loci_to_compare=[]
        for item in locus_counts.iteritems():
            if item[1] == 2:
                loci_to_compare.append(item[0])

        phase_writer = GffWriter(self.args.proj+"/annotate/phase.gff")
        phase_writer.writeMetaData('pacbio-variant-version', '1.4')

        for locus in loci_to_compare:
            seqs_to_compare=[]
            for item in self.locus_dict.iteritems():
                if item[1] == locus:
                    seqs_to_compare.append(item[0]+"|quiver")
            sequences = extract_sequence(self.args.proj+"/reseq/reseq_output_1.fasta", seqs_to_compare)
            write_fasta( sequences, tmp_fn, "w")
            ### we will align both seqs from the same locus to the MSA profile so that we can compare them
            output = self.args.proj+"/annotate/locus_"+locus+"_comparison.afa"
            MSA_fn = MSA_fn_dict[locus]
            if not os.path.isfile(output):
                try:
                    muscle_output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
                        muscle -in %s -out %s 2> /dev/null" % (tmp_fn, output),
                        executable='/bin/bash', shell=True)
                except:
                    self.log.info("MSA failed for ( %s )" % (locus) )
                    os.remove(tmp_fn)
                    continue
            sequences = extract_sequence(output, seqs_to_compare)
            if len(sequences) != 2: continue
            alignment = MSA_aligner(sequences[0].sequence, sequences[1].sequence, mode = 'string_to_string')
            qname = sequences[0].name
            tname = sequences[1].name
            q_vars = create_var_annotation(sequences[0].sequence, alignment.qvars, qname )          
            for record in q_vars:
                phase_writer.writeRecord(record)
            t_vars = create_var_annotation(sequences[1].sequence, alignment.tvars, tname )
            for record in t_vars:
                phase_writer.writeRecord(record)
        return
                
"""
    def phase_unmapped_reads( self ):
        create_directory( self.args.proj+"/unmapped" )
        num_unmapped_reads = fasta_size(fasta_fn)
        self.log.info("There are ( %s ) unmapped reads, ( %s ) percent will be passed to HGAP." % ( num_unmapped_reads, self.args.HGAP_dilute) )
        f = FastaReader(fasta_fn)
        with open(self.args.proj+"/unmapped/HGAP_input.fasta", "w") as of:
        for r in f:
            if random.random() < float(self.args.HGAP_dilute):
            print>>of, ">"+r.name
            print>>of, r.sequence
        fasta_fn = self.args.proj+"/unmapped/HGAP_input.fasta"
            
        ### these are the options that get passed to HGAP
        ### ive chosen not to make them available at the command line because theres too damn many
        ### and I dont know what most of them do anyway
        HGAP_config = """
"""
[General]
input_fofn = %s
input_type = fasta
length_cutoff = 2500
RQ_threshold = 0.75
sge_option_dm = -pe smp 4 -q secondary
sge_option_pa = -pe smp 6 -q secondary
sge_option_ca = -pe smp 4 -q secondary
sge_option_qv = -pe smp 16 -q secondary
sge_option_ck = -pe smp 1 -q secondary 
SEYMOUR_HOME = /mnt/secondary/Smrtpipe/builds/Assembly_Mainline_Nightly_Archive/build470-116466/
bestn = 48
target = pre_assembly
preassembly_num_chunk = 1 
dist_map_num_chunk = 1 
tmpdir = /tmp
big_tmpdir = /tmp
sge_queue = secondary
min_cov = 8
max_cov = 100
trim_align = 50
trim_align = 50
chunk_bestn = 5
q_nproc = 4 
""""""     % ( self.args.proj+"/unmapped/input.fofn" )     
        with open(self.args.proj+"/unmapped/HGAP.cfg", "w") as of: print>>of, HGAP_config
        with open(self.args.proj+"/unmapped/input.fofn", "w") as of: print>>of, fasta_fn
        ### this is the HGAp I hacked to allow a fasta as input
        HGAP_executable = '/home/UNIXHOME/jquinn/HLA/HGAp/HGAP_JQ.py' 
        self.log.info("HGAP envoked")
        output = check_output( "source /home/UNIXHOME/jchin/Share/HGA_env/bin/activate; cd %s/unmapped; %s %s" % (self.args.proj, HGAP_executable, self.args.proj+"/unmapped/HGAP.cfg"), executable='/bin/bash', shell=True )       
        self.log.info("HGAP exited with the following output:\n%s" % output )
        runbash("cat %s > %s" % (self.args.proj+"/unmapped/dist_map/pre_assembled_reads_*.fa", self.args.proj+"/unmapped/pre_assembled_reads.fasta") )
        self.log.info("HGAP created ( %s ) PLRs." % ( fasta_size(self.args.proj+"/unmapped/pre_assembled_reads.fasta") ) )
        trim_fasta(self.args.proj+"/unmapped/pre_assembled_reads.fasta", 
            self.args.proj+"/unmapped/pre_assembled_reads_trimmed.fasta", 
            mode="redundant", 
            pctid_threshold=float(90), 
            coverage=None, 
            aln_portion=float(0.50) )   
        self.log.info("After trimming there are ( %s ) unique PLRs" % ( fasta_size(self.args.proj+"/unmapped/pre_assembled_reads_trimmed.fasta") ) )
"""
if __name__ == '__main__':
    pipeline = HlaPipeline()
    pipeline()
