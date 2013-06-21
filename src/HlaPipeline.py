#! /usr/bin/env python
import re, os, sys, csv
import shutil
import random
import math
import logging

from subprocess import check_output
from collections import namedtuple, Counter

from pbcore.io.GffIO import GffReader, Gff3Record, GffWriter
from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord
from pbcore.io.FastqIO import FastqReader, FastqWriter

from pbhla.options import args, parse_args
from pbhla.separate_sequences import (separate_sequences,
                                      separate_listed_sequences, 
                                      separate_aligned_sequences, 
                                      write_group,
                                      write_all_groups)
from pbhla.fasta.utils import (write_fasta, 
                               fasta_size, 
                               extract_sequence, 
                               combine_fasta)
from pbhla.io.SubreadExtractor import SubreadExtractor
from pbhla.io.SamIO import SamReader
from pbhla.io.GffIO import create_annotation, create_var_annotation
from pbhla.reference.ReferenceDict import ReferenceDict
from pbhla.reference.ReferenceSelector import ReferenceSelector
from pbhla.reference.LocusReference import LocusReference
from pbhla.stats.AmpliconFinder import AmpliconFinder
from pbhla.stats.SubreadStats import SubreadStats
from pbhla.align.MultiSequenceAligner import MSA_aligner
from pbhla.phasing.clusense import Clusense
from pbhla.phasing.SummarizeClusense import combine_clusense_output
from pbhla.resequencing.Resequencer import resequence
from pbhla.annotate.summarize import (summarize_contigs,
                                      meta_summarize_contigs)
from pbhla.external.HbarTools import HbarRunner
from pbhla.external.BlasrRunner import BlasrRunner
from pbhla.external.CdHit import cd_hit_est
from pbhla.utils import *

log = logging.getLogger()

class HlaPipeline( object ):

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
                  'resequencing', 'results', 'stats']:
            sub_dir = os.path.join( args.output, d )
            create_directory( sub_dir )
            setattr(self, d, sub_dir)

    def _initialize_logging( self ):
        log_format = "%(asctime)s [%(levelname)s] %(funcName)s %(message)s"
        log_file = os.path.join( self.log_files, "HLA_Pipeline.log" )
        logging.basicConfig( level=logging.INFO, 
                             format=log_format )
        logging.basicConfig( level=logging.INFO, 
                             format=log_format, 
                             filename=log_file )
    
    def run( self ):
        # First we assemble the supplied data via HGAP / HBAR
        self.extract_subreads()
        self.create_input_fofn()
        #self.parse_fullpass_readlist()
        self.assemble_contigs()
        self.export_hbar_subreads()
        self.create_renaming_key()
        # Second we extract the subreads ourselves and map them onto
        #     the created contigs
        self.align_subreads_to_contigs()
        # Third, we align the contigs to the Human Genome ...
        self.align_contigs_to_genome()
        # ... and the locus reference sequences ...
        self.create_locus_reference()
        self.align_contigs_to_reference()
        # ... in order to separate the on-target from off-target hits
        self.find_on_target_contigs()
        self.separate_off_target_contigs()
        # Fourth, we combine the Subread and Contig dicts and use them
        #     to separate out off-target subreads
        self.separate_off_target_subreads()
        # Next we create re-align just the HLA data for summarizing
        self.trim_contigs()
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
        if args.resequence:
            self.resequence_contigs()
            self.separate_resequenced_contigs()
            self.realign_subreads_to_resequenced()
            self.align_resequenced_to_reference()
            self.separate_subreads_by_resequenced_contig()
            self.summarize_resequenced_by_locus()
        else:
            self.align_phased_to_reference()
            self.summarize_phased_by_locus()
        cleanup_directory( self.subreads )

    def extract_subreads( self ):
        log.info('Beginning the extraction of subread data')
        # Dump all valid reads from the above files
        output_path = os.path.join( self.subreads, "all_subreads" )
        self.raw_subread_file = output_path + '.fasta'
        self.fullpass_file = output_path + '_fullpass.txt'
        if os.path.isfile( self.raw_subread_file ):
            log.info('Found existing subread file "%s"' % self.raw_subread_file)
            log.info('Skipping subread extraction step\n')
            return
        log.info('Extracting subreads from input files...')
        extractor = SubreadExtractor( args.input_file, 
                                      output_path,
                                      min_length=args.min_read_length,
                                      min_score=args.min_read_score,
                                      dilution=args.dilution)
        extractor()
        check_output_file( self.raw_subread_file )
        check_output_file( self.fullpass_file )
        log.info('Finished the extraction of subread data\n')

    def create_input_fofn(self):
        log.info('Creating input fofn')
        self.input_fofn = os.path.abspath( os.path.join( self.HBAR, 'input.fofn') )
        if os.path.isfile( self.input_fofn ):
            log.info('Existing input fofn file found, skipping...\n')
            return
        with open( self.input_fofn, 'w' ) as handle:
            handle.write( os.path.abspath( self.raw_subread_file ))

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
        hbar = HbarRunner( self.input_fofn, 
                           self.HBAR,
                           min_length=args.min_read_length,
                           min_score=args.min_read_score)
        hbar()
        # Copy the contig file to a more convenient location
        shutil.copy( contig_output, self.contig_file )
        check_output_file( self.contig_file )
        log.info('Finished the assembly of subread data\n')

    def export_hbar_subreads(self):
        log.info('Exporting renamed sub-reads from HBAR')
        self.subread_file = os.path.join( self.subreads, 'renamed_subreads.fasta' )
        if os.path.isfile( self.subread_file ):
            log.info('Found existing renamed subread file "%s"' % self.subread_file)
            log.info('Skipping HBAR subread export step\n')
            return
        # Copy the renamed files locally
        fasta_folder = os.path.join( self.HBAR, '0-fasta_files' )
        for entry in os.listdir( fasta_folder ):
            if entry.endswith('_q.fa'):
                hbar_fasta = os.path.join( fasta_folder, entry )
                shutil.copy( hbar_fasta, self.subread_file )
        log.info('Finished exporting HBAR subreads\n')

    def create_renaming_key(self):
        self.renaming_key = os.path.join( self.subreads, 'renaming_key.txt' )
        if os.path.isfile( self.renaming_key ):
            log.info('Found existing renaming key "%s"' % self.renaming_key)
            log.info('Skipping key generation step\n')
            return
        # Compare the two files to make sure they're equivalent
        raw_count = fasta_size( self.raw_subread_file )
        new_count = fasta_size( self.subread_file )
        try:
            assert raw_count == new_count 
        except AssertionError:
            msg = 'The number of raw subreads (%s) does not ' % raw_count + \
                  'match the number of renamed reads (%s)' % new_count
            log.info( msg )
            raise ValueError( msg )
        # Write out the pairs of names to file
        raw_subreads = FastaReader( self.raw_subread_file )
        renamed_subreads = FastaReader( self.subread_file )
        with open( self.renaming_key, 'w') as handle:
            for raw, renamed in zip(raw_subreads, renamed_subreads):
                raw_name = raw.name.split()[0]
                new_name = renamed.name.split()[0]
                handle.write('%s\t%s\n' % (new_name, raw_name))
        self.renaming_dict = read_dict_file( self.renaming_key )
        print self.renaming_dict

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
        else:
            subread_count = fasta_size( self.subread_file )
            contig_count = fasta_size( self.contig_file )
            # Run BLASR
            log.info("Aligning {0} contigs to the {1} contigs".format(subread_count, contig_count))
            blasr_args = {'nproc': args.nproc,
                          'bestn': 1,
                          'nCandidates': 30}
            BlasrRunner( self.subread_file, self.contig_file, self.subreads_to_contigs, blasr_args )
            # Check and save the output
            check_output_file( self.subreads_to_contigs )
        log.info("Finished aligning subreads to the created contigs\n")
        self.subread_contig_dict = ReferenceDict( self.subreads_to_contigs )

    def align_contigs_to_genome(self):
        log.info('"Aligning all contigs to the Human Genome')
        self.contigs_to_genome = os.path.join(self.alignments, 'contigs_to_genome.m1')   
        if os.path.isfile( self.contigs_to_genome ):
            log.info('Found existing alignment file "{0}"'.format(self.contigs_to_genome))
            log.info("Skipping alignment step...\n")
        else:
            sa_file = args.genome + '.sa'
            contig_count = fasta_size( self.contig_file )
            # Run BLASR
            log.info("Aligning {0} contigs to the genomic reference".format(contig_count))
            blasr_args = {'nproc': args.nproc,
                          'sa': sa_file,
                          'bestn': 1,
                          'nCandidates': 10}
            BlasrRunner( self.contig_file, args.genome, self.contigs_to_genome, blasr_args )
            # Check the output and convert it to 
            check_output_file( self.contigs_to_genome )
        log.info("Finished aligning contigs to the genomic reference\n")
        self.contig_genome_dict = ReferenceDict( self.contigs_to_genome )
        self.reference_locus_dict = ReferenceDict( args.reference_file )

    def create_locus_reference(self):
        log.info("Creating locus reference fasta file")
        self.reference_seqs = os.path.join(self.references, "locus_references.fasta")
        # If no locus key and Choose_Ref is None, read the locus from the regular reference
        if os.path.isfile( self.reference_seqs ):
            log.info('Found existing locus reference sequences "{0}"'.format(self.reference_seqs))
            log.info('Skipping reference extraction step\n')
            return
        log.info("No locus reference sequences found, creating one...")
        LocusReference(args.reference_file, self.reference_seqs)
        check_output_file( self.reference_seqs )
        log.info("Finished creating locus reference\n")
    
    def align_contigs_to_reference(self):
        log.info('"Aligning all HLA contigs to the Reference sequences')
        self.contigs_to_reference = os.path.join(self.alignments, 'contigs_to_reference.m1')   
        if os.path.isfile( self.contigs_to_reference ):
            log.info('Found existing alignment file "{0}"'.format(self.contigs_to_reference))
            log.info("Skipping alignment step...\n")
        else:   
            contig_count = fasta_size( self.contig_file )
            reference_count = fasta_size( self.reference_seqs )
            # Run BLASR
            log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
            blasr_args = {'nproc': args.nproc,
                          'bestn': 1,
                          'nCandidates': reference_count}
            BlasrRunner( self.contig_file, self.reference_seqs, self.contigs_to_reference, blasr_args )
            # Check the output and convert it into a dictionary
            check_output_file( self.contigs_to_reference )
        log.info("Finished aligning contigs to the HLA reference set\n")
        self.contig_reference_dict = ReferenceDict( self.contigs_to_reference )
        self.contig_locus_dict = cross_ref_dict( self.contig_reference_dict, self.reference_locus_dict )

    def find_on_target_contigs(self):
        log.info('Identifying on-target contigs')
        self.on_target_contigs = os.path.join( self.references, 'hla_contig_ids.txt' )
        on_target_contig_ids = [c for c in self.contig_locus_dict
                                   if self.contig_genome_dict[c] == 'chr6']
        with open(self.on_target_contigs, 'w') as handle:
            for contig_id in on_target_contig_ids:
                handle.write( contig_id + '\n' )
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
        separated_seqs = separate_listed_sequences( self.contig_file, 
                                                    self.on_target_contigs )
        write_group(separated_seqs, 'selected', self.hla_contigs)
        write_group(separated_seqs, 'not_selected', self.non_hla_contigs)
        log.info('Finished separating off-target contigs\n')
    
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
        separated_seqs = separate_aligned_sequences( self.subread_file, 
                                                     self.subread_contig_dict,
                                                     self.on_target_contigs)
        write_group(separated_seqs, 'selected', self.hla_subreads)
        write_group(separated_seqs, 'not_selected', self.non_hla_subreads)
        log.info('Finished separating off-target subreads\n')

    def trim_contigs(self):
        log.info('Trimming contigs to their aligned region')
        self.trimmed_contigs = os.path.join(self.references, 'trimmed_contigs.fasta')
        if os.path.isfile( self.trimmed_contigs ):
            log.info('Found existing Trimmed contig file "%s"' % self.trimmed_contigs )
            log.info('Skipping contig trimming step...\n')
            return
        trims = {}
        with open(self.contigs_to_reference) as handle:
            for record in map(BlasrM1._make, csv.reader(handle, delimiter=' ')):
                start = max(int(record.qstart)-100, 0)
                end = min(int(record.qend)+100, int(record.qlength))
                trims[record.qname] = (start, end)
        with FastaWriter(self.trimmed_contigs) as writer:
            for record in FastaReader(self.hla_contigs):
                name = record.name.split()[0]
                start, end = trims[name]
                record = FastaRecord( name, record.sequence[start:end] )
                writer.writeRecord( record )
        log.info('Finished trimming contig sequences\n')

    def remove_redundant_contigs(self):
        log.info("Removing redundant contigs from the HLA contig file")
        self.selected_contigs = os.path.join(self.references, 'selected_contigs.fasta')
        if os.path.isfile( self.selected_contigs ):
            log.info('Found existing non-redundant contig file "{0}"'.format(self.selected_contigs))
            log.info("Skipping redundant contig filtering step...\n")
            return
        log.info('File of selected contigs not found, initializing Cd-Hit-Est')
        cd_hit_est( self.trimmed_contigs, self.selected_contigs )
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
        blasr_args = {'nproc': args.nproc,
                      'bestn': 1,
                      'nCandidates': ref_count}
        BlasrRunner(self.hla_subreads, self.selected_contigs, self.hla_alignment, blasr_args)
        check_output_file( self.hla_alignment )
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
        contig_fofn = os.path.join( self.references, "Contig_Files.txt" )
        if os.path.isfile( contig_fofn ):
            log.info('Found existing Contig File List "{0}"'.format(contig_fofn))
            log.info("Skipping contig separation step...\n")
            self.contig_files = read_list_file( contig_fofn )
            return
        separated_seqs = separate_sequences( self.selected_contigs )
        contig_prefix = os.path.join(self.references, 'Contig')
        self.contig_files = write_all_groups( separated_seqs, contig_prefix )
        write_list_file( self.contig_files, contig_fofn)
        log.info('Finished separating contigs')

    def create_subread_selected_contig_dict(self):
        log.info("Converting the Subreads-to-Contigs alignment to a Dictionary")
        self.subread_selected_contig_dict = ReferenceDict( self.hla_alignment )
        log.info("Finished convereting the data to a Dictionary\n")

    def separate_subreads_by_contig(self):
        log.info("Separating subreads by aligned contig")
        subread_fofn = os.path.join( self.subreads, "Subread_Files.txt" )
        if os.path.isfile( subread_fofn ):
            log.info('Found existing Subread File List "{0}"'.format(subread_fofn))
            log.info("Skipping subread separation step...\n")
            self.subread_files = read_list_file( subread_fofn )
            return
        separated_seqs = separate_sequences( self.hla_subreads, 
                                             self.subread_selected_contig_dict )
        subread_prefix = os.path.join(self.subreads, 'Subreads')
        self.subread_files = write_all_groups( separated_seqs, subread_prefix )
        self.subread_files = [fn for fn in self.subread_files if not fn.endswith('Unmapped.fasta')]
        write_list_file( self.subread_files, subread_fofn )
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
                Clusense(subread_fn, contig_fn, output_folder, threshold=0.1)
            cleanup_directory( output_folder )
        log.info('Finished Phasing subreads with Clusense')

    def summarize_clusense(self):
        log.info("Summarizing the output from Clusense")
        output_folder = os.path.join(self.phasing_results, 'Clusense_Results')
        self.clusense_cns_fofn, self.clusense_read_fofn = combine_clusense_output(self.phasing, output_folder)
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

    def resequence_contigs(self):
        log.info("Resequencing phased contigs")
        self.resequenced_contigs = os.path.join(self.references, 'resequenced_contigs.fasta')
        if os.path.exists( self.resequenced_contigs ):
            log.info('Existing Quiver consensus output found, skipping...')
            return
        consensus_file = resequence( args.smrt_analysis, args.input_file, self.phased_contigs, self.resequencing)
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
        else:
            contig_count = fasta_size( self.resequenced_hla )
            reference_count = fasta_size( self.reference_seqs )
            # Run BLASR
            log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
            blasr_args = {'nproc': args.nproc,
                          'bestn': 1,
                          'nCandidates': reference_count}
            BlasrRunner( self.resequenced_hla, self.reference_seqs, self.resequenced_to_reference, blasr_args )
            # Check and save the output
            check_output_file( self.resequenced_to_reference )
        self.resequenced_reference_dict = ReferenceDict( self.resequenced_to_reference )
        self.resequenced_locus_dict = cross_ref_dict( self.resequenced_reference_dict, self.reference_locus_dict )
        log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def realign_subreads_to_resequenced(self):
        log.info("Re-aligning HLA subreads to resequenced references")
        self.resequenced_alignment = os.path.join( self.alignments, "hla_subreads_to_resequenced.sam" )
        if os.path.isfile( self.resequenced_alignment ):
            log.info('Found existing SAM file "{0}"'.format(self.resequenced_alignment))
            log.info("Skipping realignment step...\n")
        else:
            query_count = fasta_size( self.hla_subreads )
            ref_count = fasta_size( self.resequenced_hla )
            log.info("Aligning {0} subreads to {1} reference sequences".format(query_count, ref_count))
            blasr_args = {'nproc': args.nproc,
                          'bestn': 1,
                          'nCandidates': ref_count}
            BlasrRunner( self.hla_subreads, 
                         self.resequenced_hla, 
                         self.resequenced_alignment, 
                         blasr_args )
            check_output_file( self.resequenced_alignment )
        self.subread_resequenced_contig_dict = ReferenceDict( self.resequenced_alignment )
        log.info("Finished realigning HLA subreads to resequenced contigs\n")

    def separate_subreads_by_resequenced_contig(self):
        log.info("Separating subreads by resequenced contig")
        self.subread_fofn = os.path.join( self.subreads, "resequenced_subread_files.txt" )
        if os.path.isfile( self.subread_fofn ):
            log.info('Found existing Subread File List "{0}"'.format(self.subread_fofn))
            log.info("Skipping subread separation step...\n")
            self.subread_files = read_list_file( self.subread_fofn )
            return
        separator = SequenceSeparator( self.hla_subreads, 
                                       reference_dict=self.subread_resequenced_contig_dict )
        subread_prefix = os.path.join(self.subreads, 'Resequenced')
        self.subread_files = separator.write_all( subread_prefix )
        self.subread_files = [fn for fn in self.subread_files if not fn.endswith('Unmapped.fasta')]
        write_list_file( self.subread_files, self.subread_fofn )
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
        else:
            contig_count = fasta_size( self.phased_hla )
            reference_count = fasta_size( self.reference_seqs )
            # Run BLASR
            log.info("Aligning {0} contigs to {1} reference sequences".format(contig_count, reference_count))
            blasr_args = {'nproc': args.nproc,
                          'bestn': 1,
                          'nCandidates': reference_count}
            BlasrRunner( self.phased_hla, self.reference_seqs, self.phased_to_reference, blasr_args )
            # Check and save the output
            check_output_file( self.phased_to_reference )
        self.phased_reference_dict = ReferenceDict( self.phased_to_reference )
        self.phased_locus_dict = cross_ref_dict( self.phased_reference_dict, self.reference_locus_dict )
        log.info("Finished aligning resequenced contigs to the HLA reference set\n")

    def summarize_phased_by_locus(self):
        log.info("Picking contigs from the resequencing and realignment results")
        clusense_results = os.path.join( self.results, 'Clusense' )
        create_directory( clusense_results )
        summaries = summarize_contigs( self.phased_hla, 
                                       self.clusense_read_fofn, 
                                       self.phased_locus_dict, 
                                       self.phased_to_reference,
                                       clusense_results )
        meta_summarize_contigs( summaries, clusense_results )
        log.info("Finished picking contigs")

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
    HlaPipeline().run()
