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
from pbhla.io.SubreadExtractor import SubreadExtractor
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
from pbhla.annotate.ContigPicker import ContigPicker
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
MIN_LENGTH = 2000
NUM_PROC = 4

log = logging.getLogger()

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
        add("--smrt_path", metavar="PATH", default=SMRT_ANALYSIS,
            help="Path to the setup script for the local SMRT Analysis installation")
        add("--resequence", action="store_true", help="Use quiver to resequence the")
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
                  'resequencing', 'results', 'stats']:
            sub_dir = os.path.join( self.output, d )
            create_directory( sub_dir )
            setattr(self, d, sub_dir)

    def initialize_logging( self ):
        time_format = "%I:%M:%S"
        log_format = "%(asctime)s %(levelname)s %(filename)s " + \
                     "%(funcName)s %(lineno)d %(message)s"
        out_format = "%(asctime)s %(filename)s %(funcName)s %(message)s"
        self.log = logging.getLogger()
        log.setLevel( logging.INFO )
        # Set-up one logger for STDOUT
        h1 = logging.StreamHandler( stream=sys.stdout )
        h1.setLevel( logging.INFO )
        f1 = logging.Formatter( fmt=out_format, datefmt=time_format )
        h1.setFormatter( f1 )
        log.addHandler( h1 )
        # Setup a second logger to log to file
        log_file = os.path.join( self.log_files, "HLA_Pipeline.log" )
        h2 = logging.FileHandler( log_file )
        f2 = logging.Formatter( fmt=log_format, datefmt=time_format )
        h2.setFormatter( f2 )
        h2.setLevel( logging.INFO )
        log.addHandler( h2 )
    
    def validate_settings( self ):
        # Report the settings with which the pipeline was invoked
        log.info( "HLA Pipeline invoked:\n\t{0}\n".format(" ".join(sys.argv)))
        self.input_file = os.path.abspath( self.input_file )
        # Check dilution factors
        if self.dilution <= 0.0 or self.dilution > 1.0:
            msg = "Dilute factor must be between 0 and 1"
            log.info( msg )
            raise ValueError( msg )
        if self.resequence:
            self.smrt_analysis = SmrtAnalysisRunner( self.smrt_path, self.log_files, self.nproc )
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

    def getVersion(self):
        return __version__

    def __call__( self ):
        self.create_input_fofn()
        # First we assemble the supplied data via HGAP / HBAR
        self.assemble_contigs()
        # Second we extract the subreads ourselves and map them onto
        #     the created contigs
        self.extract_subreads()
        self.parse_fullpass_readlist()
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
        self.find_chimeric_contigs()
        self.separate_chimeric_contigs()
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
        if self.resequence:
            self.resequence_contigs()
            self.separate_resequenced_contigs()
            self.realign_subreads_to_resequenced()
            self.align_resequenced_to_reference()
            self.create_resequenced_dict()
            self.create_resequenced_locus_dict()
            self.create_subread_resequenced_contig_dict()
            self.separate_subreads_by_resequenced_contig()
            self.summarize_resequenced_by_locus()
        else:
            self.align_phased_to_reference()
            self.create_phased_dict()
            self.create_phased_locus_dict()
            self.summarize_phased_by_locus()
        self.cleanup_subreads()

    def check_output_file(self, filepath):
        try:
            assert os.path.isfile( filepath )
        except:
            msg = 'Expected output file not found! "{0}"'.format(filepath)
            log.info( msg )
            raise IOError( msg )

    def create_input_fofn(self):
        log.info('Creating input fofn')
        self.input_fofn = os.path.abspath( os.path.join( self.output, 'input.fofn') )
        if os.path.isfile( self.input_fofn ):
            log.info('Existing input fofn file found, skipping...\n')
            return
        # Identify Bas.H5 and Bax.H5 files
        if self.input_file.endswith('.bax.h5') or self.input_file.endswith('.bas.h5'):
            bash_files = [ self.input_file ]
        if self.input_file.endswith('.fofn'):
            with open( self.input_file ) as handle:
                bash_files = [ line.strip() for line in handle ]
        # Convert the Bas.H5 files to Bax.H5 files
        bax_files = []
        for filename in bash_files:
            if filename.endswith('.bax.h5'):
                bax_files.append( filename )
            elif filename.endswith('.bas.h5'):
                root = '.'.join( filename.split('.')[:-2] )
                for i in range(1,4):
                    bax_name = '{0}.{1}.bax.h5'.format(root, i)
                    bax_files.append( bax_name )
        # Write the Bax.H5 files out
        with open(self.input_fofn, 'w') as handle:
            for filename in bax_files:
                print >> handle, os.path.abspath( filename )
        log.info('Finished creating input fofn')

    def assemble_contigs( self ):
        log.info('Beginning the extraction of subiread data')
        contig_output = os.path.join( self.HBAR,
                                      "3-CA",
                                      "9-terminator",
                                      "asm.utg.fasta" )
        self.contig_file = os.path.join( self.references,
                                         "all_contigs.fasta")
        if os.path.isfile( self.contig_file ):
            log.info('Found existing contig file "{0}"'.format(self.contig_file))
            log.info('Skipping HBAR assembly step\n')
            return
        # Run HGAP
        log.info('No contig file found, initializing new HbarRunner')
        hbar = HbarRunner( self.input_fofn, self.HBAR )
        hbar()
        # Copy the contig file to a more convenient location
        shutil.copy( contig_output, self.contig_file )
        self.check_output_file( self.contig_file )
        log.info('Finished the assembly of subread data\n')

    def extract_subreads( self ):
        log.info('Beginning the extraction of subread data')
        # Dump all valid reads from the above files
        output_path = os.path.join( self.subreads, "all_subreads" )
        self.subread_file = output_path + '.fasta'
        self.fullpass_file = output_path + '_fullpass.txt'
        if os.path.isfile( self.subread_file ):
            log.info('Found existing subread file "{0}"'.format(self.subread_file))
            log.info('Skipping subread extraction step\n')
            return
        log.info('Extracting subreads from input files...')
        extractor = SubreadExtractor( self.input_fofn, 
                                      output_path,
                                      min_length=self.min_read_length,
                                      min_score=self.min_read_score,
                                      dilution=self.dilution)
        extractor()
        self.check_output_file( self.subread_file )
        self.check_output_file( self.fullpass_file )
        log.info('Finished the extraction of subread data\n')

    def parse_fullpass_readlist(self):
        log.info('Reading the list of fullpass reads...')
        fullpass_reads = []
        with open( self.fullpass_file, 'r') as handle:
            for line in handle:
                fullpass_reads.append( line.strip() )
        self.fullpass_reads = set(fullpass_reads)
        log.info('Finished parsing the fullpass readlist\n')

    def align_subreads_to_contigs(self):
        log.info('"Aligning all subreads to the HBAR contigs')
        self.subreads_to_contigs = os.path.join(self.alignments, 'subreads_to_contigs.m1')   
        if os.path.isfile( self.subreads_to_contigs ):
            log.info('Found existing alignment file "{0}"'.format(self.subreads_to_contigs))
            log.info("Skipping alignment step...\n")
            return
        subread_count = fasta_size( self.subread_file )
        contig_count = fasta_size( self.contig_file )
        # Run BLASR
        log.info("Aligning {0} contigs to the {1} contigs".format(subread_count, contig_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': 30}
        BlasrRunner( self.subread_file, self.contig_file, self.subreads_to_contigs, blasr_args )
        # Check and save the output
        self.check_output_file( self.subreads_to_contigs )
        log.info("Finished aligning subreads to the created contigs\n")

    def create_subread_contig_dict(self):
        log.info("Converting the Subreads-to-Contigs alignment to a Dictionary")
        self.subread_contig_dict = ReferenceDict( self.subreads_to_contigs )
        log.info("Finished convereting the data to a Dictionary\n")

    def align_contigs_to_genome(self):
        log.info('"Aligning all contigs to the Human Genome')
        self.contigs_to_genome = os.path.join(self.alignments, 'contigs_to_genome.m1')   
        if os.path.isfile( self.contigs_to_genome ):
            log.info('Found existing alignment file "{0}"'.format(self.contigs_to_genome))
            log.info("Skipping alignment step...\n")
            return
        sa_file = self.genome + '.sa'
        contig_count = fasta_size( self.contig_file )
        # Run BLASR
        log.info("Aligning {0} contigs to the genomic reference".format(contig_count))
        blasr_args = {'nproc': self.nproc,
                      'sa': sa_file,
                      'bestn': 1,
                      'nCandidates': 10}
        BlasrRunner( self.contig_file, self.genome, self.contigs_to_genome, blasr_args )
        # Check and save the output
        self.check_output_file( self.contigs_to_genome )
        log.info("Finished aligning contigs to the genomic reference\n")

    def create_contig_genome_dict(self):
        log.info("Converting the Contigs-to-Genome alignment to a Dictionary")
        self.contig_genome_dict = ReferenceDict( self.contigs_to_genome )
        log.info("Finished convereting the data to a Dictionary\n")

    def create_ref_locus_dict(self):
        log.info('"Creating reference dictonary from "{0}"'.format(self.reference_file))
        self.reference_locus_dict = ReferenceDict( self.reference_file )
        log.info("Finished creating locus reference dict...\n")

    def create_locus_reference(self):
        log.info("Creating locus reference fasta file")
        self.reference_seqs = os.path.join(self.references, "locus_references.fasta")
        # If no locus key and Choose_Ref is None, read the locus from the regular reference
        if os.path.isfile( self.reference_seqs ):
            log.info('Found existing locus reference sequences "{0}"'.format(self.reference_seqs))
            log.info('Skipping reference extraction step\n')
            return
        log.info("No locus reference sequences found, creating one...")
        LocusReference(self.reference_file, self.reference_seqs)
        self.check_output_file( self.reference_seqs )
        log.info("Finished creating locus reference\n")
    
    def align_contigs_to_reference(self):
        log.info('"Aligning all HLA contigs to the Reference sequences')
        self.contigs_to_reference = os.path.join(self.alignments, 'contigs_to_reference.m1')   
        if os.path.isfile( self.contigs_to_reference ):
            log.info('Found existing alignment file "{0}"'.format(self.contigs_to_reference))
            log.info("Skipping alignment step...\n")
            return
        contig_count = fasta_size( self.contig_file )
        reference_count = fasta_size( self.reference_seqs )
        # Run BLASR
        log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': reference_count}
        BlasrRunner( self.contig_file, self.reference_seqs, self.contigs_to_reference, blasr_args )
        # Check and save the output
        self.check_output_file( self.contigs_to_reference )
        log.info("Finished aligning contigs to the HLA reference set\n")

    def create_contig_reference_dict(self):
        log.info("Converting the Contigs-to-Reference alignment to a Dictionary")
        self.contig_reference_dict = ReferenceDict( self.contigs_to_reference )
        log.info("Finished convereting the data to a Dictionary\n")

    def create_contig_locus_dict(self):
        log.info("Converting the Contigs-to-Reference dict and Reference-to-Locus")
        log.info("    dict into a combined Contig-to-Locus Dictionary")
        self.contig_locus_dict = cross_ref_dict( self.contig_reference_dict, self.reference_locus_dict )
        log.info("Finished convereting the data to a Dictionary\n")

    def find_on_target_contigs(self):
        log.info('Identifying on-target contigs')
        self.on_target_contig_ids = [c for c in self.contig_locus_dict
                                       if self.contig_genome_dict[c] == 'chr6']
        log.info('Finished identifying on-target contigs\n')

    def separate_off_target_contigs(self):
        log.info('Separating off-target contigs by genomic alignment')
        self.hla_contigs = os.path.join(self.references, 'hla_contigs.fasta')
        self.non_hla_contigs = os.path.join(self.references, 'non_hla_contigs.fasta')
        if os.path.isfile( self.hla_contigs ) and os.path.isfile( self.non_hla_contigs ):
            log.info('Found existing sequence file "{0}"'.format(self.hla_contigs))
            log.info('Found existing sequence file "{0}"'.format(self.non_hla_contigs))
            log.info("Skipping separation step...\n")
            return
        log.info('No separated contig files found, initializing separator')
        separator = SequenceSeparator( self.contig_file, selected=self.on_target_contig_ids )
        separator.write('selected', self.hla_contigs)
        separator.write('not_selected', self.non_hla_contigs)
        log.info('Finished separating off-target contigs\n')
    
    def find_on_target_subreads(self):
        log.info('Identifying on-target subreads')
        self.on_target_subreads = [s for s in self.subread_contig_dict
                                   if self.subread_contig_dict[s] in self.on_target_contig_ids]
        log.info('Finished identifying on-target subreads\n')

    def separate_off_target_subreads(self):
        log.info('Separating off-target subreads by contig alignment')
        self.hla_subreads = os.path.join(self.subreads, 'hla_subreads.fasta')
        self.non_hla_subreads = os.path.join(self.subreads, 'non_hla_subreads.fasta')
        if os.path.isfile( self.hla_subreads ) and os.path.isfile( self.non_hla_subreads ):
            log.info('Found existing sequence file "{0}"'.format(self.hla_subreads))
            log.info('Found existing sequence file "{0}"'.format(self.non_hla_subreads))
            log.info("Skipping separation step...\n")
            return
        log.info('No separated contig files found, initializing separator')
        separator = SequenceSeparator( self.subread_file, selected=self.on_target_subreads )
        separator.write('selected', self.hla_subreads)
        separator.write('not_selected', self.non_hla_subreads)
        log.info('Finished separating off-target subreads\n')

    def find_chimeric_contigs(self):
        self.chimeric_ids = []
        log.info('Identifying chimeric contigs')
        for record in FastaReader( self.hla_contigs ):
            name = record.name.split()[0]
            if self.contig_locus_dict[name] in ['A', 'B', 'C']:
                if len(record.sequence) > 3600:
                    self.chimeric_ids.append( name )
                if len(record.sequence) < 3000:
                    self.chimeric_ids.append( name )
        log.info('Finished identifying chimeric contigs\n')

    def separate_chimeric_contigs(self):
        log.info('Separating off-target contigs by genomic alignment')
        self.chimeric_hla = os.path.join(self.references, 'chimeric_hla_contigs.fasta')
        self.non_chimeric_hla = os.path.join(self.references, 'non_chimeric_hla_contigs.fasta')
        if os.path.isfile( self.chimeric_hla ) and os.path.isfile( self.non_chimeric_hla ):
            log.info('Found existing sequence file "{0}"'.format(self.chimeric_hla))
            log.info('Found existing sequence file "{0}"'.format(self.non_chimeric_hla))
            log.info("Skipping separation step...\n")
            return
        log.info('No separated contig files found, initializing separator')
        separator = SequenceSeparator( self.hla_contigs, selected=self.chimeric_ids )
        separator.write('selected', self.chimeric_hla)
        separator.write('not_selected', self.non_chimeric_hla)
        log.info('Finished separating off-target contigs\n')

    def remove_redundant_contigs(self):
        log.info("Removing redundant contigs from the HLA contig file")
        self.selected_contigs = os.path.join(self.references, 'selected_contigs.fasta')
        if os.path.isfile( self.selected_contigs ):
            log.info('Found existing non-redundant contig file "{0}"'.format(self.selected_contigs))
            log.info("Skipping redundant contig filtering step...\n")
            return
        log.info('File of selected contigs not found, initializing Cd-Hit-Est')
        cd_hit_est( self.non_chimeric_hla, self.selected_contigs )
        log.info('Finished selecting non-redundant contigs\n')

    def realign_hla_subreads(self):
        log.info("Re-aligning HLA subreads to selected references")
        self.hla_alignment = os.path.join( self.alignments, "hla_subreads_to_contigs.sam" )
        if os.path.isfile( self.hla_alignment ):
            log.info('Found existing SAM file "{0}"'.format(self.hla_alignment))
            log.info("Skipping realignment step...\n")
            return
        query_count = fasta_size( self.hla_subreads )
        ref_count = fasta_size( self.selected_contigs )
        log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': ref_count}
        BlasrRunner(self.hla_subreads, self.selected_contigs, self.hla_alignment, blasr_args)
        self.check_output_file( self.hla_alignment )
        log.info("Finished realigning HLA subreads to HLA contigs\n")

    def summarize_aligned_hla_subreads(self):
        log.info("Summarizing coverage of the HLA contigs by the HLA subreads")
        self.locus_stats = os.path.join(self.stats, "locus_statistics.csv")
        self.reference_stats = os.path.join(self.stats, "reference_statistics.csv")
        ### will go through sam file, sort out the raw reads, and tabulate statistics while we do it
        if os.path.isfile( self.locus_stats ) and \
           os.path.isfile( self.reference_stats ):
            log.info('Found existing Locus Statistics at "{0}"'.format(self.locus_stats))
            log.info('Found existing Reference Statistics at "{0}"'.format(self.reference_stats))
            log.info("Skipping alignment summary step...\n")
            return
        stats = SubreadStats( self.hla_contigs, self.contig_locus_dict )
        stats.add_sam_file( self.hla_alignment )
        ### finally write out the subread statistics
        stats.write( self.locus_stats, 'locus' )
        stats.write( self.reference_stats, 'reference' )
        log.info("Finished summarizing coverage of the HLA contigs")

    def separate_contigs(self):
        log.info("Separating remaining contigs into individual files")
        contig_fofn = os.path.join( self.subreads, "contig_files.txt" )
        if os.path.isfile( contig_fofn ):
            log.info('Found existing Contig File List "{0}"'.format(contig_fofn))
            log.info("Skipping contig separation step...\n")
            self.contig_files = read_fofn( contig_fofn )
            return
        separator = SequenceSeparator( self.selected_contigs )
        contig_prefix = os.path.join(self.references, 'Contig')
        self.contig_files = separator.write_all( contig_prefix )
        write_fofn( self.contig_files, contig_fofn)
        log.info('Finished separating contigs')

    def create_subread_selected_contig_dict(self):
        log.info("Converting the Subreads-to-Contigs alignment to a Dictionary")
        self.subread_selected_contig_dict = ReferenceDict( self.hla_alignment )
        log.info("Finished convereting the data to a Dictionary\n")

    def separate_subreads_by_contig(self):
        log.info("Separating subreads by aligned contig")
        subread_fofn = os.path.join( self.subreads, "subread_files.txt" )
        if os.path.isfile( subread_fofn ):
            log.info('Found existing Subread File List "{0}"'.format(subread_fofn))
            log.info("Skipping subread separation step...\n")
            self.subread_files = read_fofn( subread_fofn )
            return
        separator = SequenceSeparator( self.hla_subreads, 
                                       reference_dict=self.subread_selected_contig_dict )
        subread_prefix = os.path.join(self.subreads, 'Subreads')
        self.subread_files = separator.write_all( subread_prefix )
        self.subread_files = [fn for fn in self.subread_files if not fn.endswith('Unmapped.fasta')]
        write_fofn( self.subread_files, subread_fofn )
        log.info('Finished separating subreads by contig')

    def phase_reads_with_clusense(self):
        log.info("Phasing subreads with Clusense")
        for subread_fn, contig_fn in zip(self.subread_files, self.contig_files):
            folder_name = os.path.basename(contig_fn).split('.')[0]
            output_folder = os.path.join(self.phasing, folder_name)
            if os.path.exists( output_folder ):
                log.info('"Clusense output detected for "{0}", skipping...'.format(folder_name))
            else:
                log.info("Phasing subreads for {0}".format(folder_name))
                Clusense(subread_fn, contig_fn, output_folder)
        log.info('Finished Phasing subreads with Clusense')

    def summarize_clusense(self):
        log.info("Summarizing the output from Clusense")
        output_folder = os.path.join(self.phasing_results, 'Clusense_Results')
        self.clusense_cns_fofn, self.clusense_read_fofn = combine_clusense_output(self.phasing, output_folder)
        log.info('Finished Phasing subreads with Clusense')

    def output_phased_contigs(self):
        log.info("Creating a combined reference of non-HLA and phased-HLA contigs")
        self.phased_contigs = os.path.join(self.references, 'phased_contigs.fasta')
        clusense_cns_files = read_fofn( self.clusense_cns_fofn )
        fasta_files = list(clusense_cns_files) + [self.non_hla_contigs]
        combine_fasta(fasta_files, self.phased_contigs)
        self.phased_hla = os.path.join(self.references, 'phased_hla.fasta')
        combine_fasta(clusense_cns_files, self.phased_hla)
        log.info('Finished creating combined reference')

    def resequence_contigs(self):
        log.info("Resequencing phased contigs")
        self.resequenced_contigs = os.path.join(self.references, 'resequenced_contigs.fasta')
        if os.path.exists( self.resequenced_contigs ):
            log.info('Existing Quiver consensus output found, skipping...')
            return
        consensus_file = resequence( self.smrt_analysis, self.input_file, self.phased_contigs, self.resequencing)
        shutil.copy( consensus_file, self.resequenced_contigs )
        log.info('Finished the contig resequencing step')

    def separate_resequenced_contigs(self):
        log.info("Resequencing phased contigs")
        self.resequenced_hla = os.path.join(self.references, 'resequenced_hla.fasta')
        self.resequenced_non_hla = os.path.join(self.references, 'resequenced_non_hla.fasta')
        if os.path.isfile( self.resequenced_hla ):
            log.info('Found existing file of resequenced HLA contigs')
            log.info('Skipping resequenced contig separation step...\n')
            return
        with FastaWriter(self.resequenced_hla) as hla:
            with FastaWriter(self.resequenced_non_hla) as non_hla:
                for record in FastaReader( self.resequenced_contigs ):
                    if record.name.startswith('Contig_'):
                        hla.writeRecord( record )
                    else:
                        non_hla.writeRecord( record )
        log.info('Finished separating the resequenced contigs')

    def align_resequenced_to_reference(self):
        log.info('"Aligning all resequenced HLA contigs to the Reference sequences')
        self.resequenced_to_reference = os.path.join(self.alignments, 'resequenced_to_reference.m1')   
        if os.path.isfile( self.resequenced_to_reference ):
            log.info('Found existing alignment file "{0}"'.format(self.resequenced_to_reference))
            log.info("Skipping alignment step...\n")
            return
        contig_count = fasta_size( self.resequenced_hla )
        reference_count = fasta_size( self.reference_seqs )
        # Run BLASR
        log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': reference_count}
        BlasrRunner( self.resequenced_hla, self.reference_seqs, self.resequenced_to_reference, blasr_args )
        # Check and save the output
        self.check_output_file( self.resequenced_to_reference )
        log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def create_resequenced_dict(self):
        log.info("Converting the Resequenced-to-Reference alignment to a Dictionary")
        self.resequenced_reference_dict = ReferenceDict( self.resequenced_to_reference )
        log.info("Finished convereting the data to a Dictionary\n")

    def create_resequenced_locus_dict(self):
        log.info("Converting the Contigs-to-Reference dict and Reference-to-Locus")
        log.info("    dict into a combined Contig-to-Locus Dictionary")
        self.resequenced_locus_dict = cross_ref_dict( self.resequenced_reference_dict, self.reference_locus_dict )
        log.info("Finished convereting the data to a Dictionary\n")

    def realign_subreads_to_resequenced(self):
        log.info("Re-aligning HLA subreads to resequenced references")
        self.resequenced_alignment = os.path.join( self.alignments, "hla_subreads_to_resequenced.sam" )
        if os.path.isfile( self.resequenced_alignment ):
            log.info('Found existing SAM file "{0}"'.format(self.resequenced_alignment))
            log.info("Skipping realignment step...\n")
            return
        query_count = fasta_size( self.hla_subreads )
        ref_count = fasta_size( self.resequenced_hla )
        log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': ref_count}
        BlasrRunner(self.hla_subreads, self.resequenced_hla, self.resequenced_alignment, blasr_args)
        self.check_output_file( self.resequenced_alignment )
        log.info("Finished realigning HLA subreads to resequenced contigs\n")

    def create_subread_resequenced_contig_dict(self):
        log.info("Converting the Subreads-to-Resequenced-Contigs alignment to a Dictionary")
        self.subread_resequenced_contig_dict = ReferenceDict( self.resequenced_alignment )
        log.info("Finished convereting the data to a Dictionary\n")

    def separate_subreads_by_resequenced_contig(self):
        log.info("Separating subreads by resequenced contig")
        self.subread_fofn = os.path.join( self.subreads, "resequenced_subread_files.txt" )
        if os.path.isfile( self.subread_fofn ):
            log.info('Found existing Subread File List "{0}"'.format(self.subread_fofn))
            log.info("Skipping subread separation step...\n")
            self.subread_files = read_fofn( self.subread_fofn )
            return
        separator = SequenceSeparator( self.hla_subreads, 
                                       reference_dict=self.subread_resequenced_contig_dict )
        subread_prefix = os.path.join(self.subreads, 'Resequenced')
        self.subread_files = separator.write_all( subread_prefix )
        self.subread_files = [fn for fn in self.subread_files if not fn.endswith('Unmapped.fasta')]
        write_fofn( self.subread_files, self.subread_fofn )
        log.info('Finished separating subreads by contig')

    def summarize_resequenced_by_locus(self):
        log.info("Picking contigs from the resequencing and realignment results")
        quiver_results = os.path.join( self.results, 'Quiver' )
        create_directory( quiver_results )
        ContigPicker( self.resequenced_hla, self.subread_fofn, self.resequenced_locus_dict, quiver_results )
        log.info("Finished picking contigs")

    def align_phased_to_reference(self):
        log.info('"Aligning all phased HLA contigs to the Reference sequences')
        self.phased_to_reference = os.path.join(self.alignments, 'phased_to_reference.m1')   
        if os.path.isfile( self.phased_to_reference ):
            log.info('Found existing alignment file "{0}"'.format(self.phased_to_reference))
            log.info("Skipping alignment step...\n")
            return
        contig_count = fasta_size( self.phased_hla )
        reference_count = fasta_size( self.reference_seqs )
        # Run BLASR
        log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
        blasr_args = {'nproc': self.nproc,
                      'bestn': 1,
                      'nCandidates': reference_count}
        BlasrRunner( self.phased_hla, self.reference_seqs, self.phased_to_reference, blasr_args )
        # Check and save the output
        self.check_output_file( self.phased_to_reference )
        log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def create_phased_dict(self):
        log.info("Converting the Resequenced-to-Reference alignment to a Dictionary")
        self.phased_reference_dict = ReferenceDict( self.phased_to_reference )
        log.info("Finished convereting the data to a Dictionary\n")

    def create_phased_locus_dict(self):
        log.info("Converting the Contigs-to-Reference dict and Reference-to-Locus")
        log.info("    dict into a combined Contig-to-Locus Dictionary")
        self.phased_locus_dict = cross_ref_dict( self.phased_reference_dict, self.reference_locus_dict )
        log.info("Finished convereting the data to a Dictionary\n")

    def summarize_phased_by_locus(self):
        log.info("Picking contigs from the resequencing and realignment results")
        clusense_results = os.path.join( self.results, 'Clusense' )
        create_directory( clusense_results )
        ContigPicker( self.phased_hla, self.clusense_read_fofn, self.phased_locus_dict, clusense_results )
        log.info("Finished picking contigs")

    def cleanup_subreads(self):
        log.info("Cleaning up unused Subread files")
        for entry in os.listdir( self.subreads ):
            if entry.endswith('aln') or entry.endswith('aln_unsorted'):
                os.remove( os.path.join( self.subreads, entry) )
        log.info("Finished clean-up process")

    ### Other functions, not currently used ###
    def find_amplicons(self):
        log.info("Finding amplicon locations within each reference")
        self.amplicon_summary = self.stats_dir + "/amplicon_summary.csv"
        if os.path.isfile( self.amplicon_summary ):
            log.info('Found existing Amplicon Summary "{0}" )'.format(self.amplicon_summary))
            log.info("Skipping amplicon finding step...\n")
            return
        AmpliconFinder(self.sam_file, self.amplicon_summary, self.locus_dict)
        log.info("Finished finding amplicon locations\n")

if __name__ == '__main__':
    pipeline = HlaPipeline()
    pipeline()
