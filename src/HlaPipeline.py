#!/home/UNIXHOME/jquinn/HGAP_env/bin/python
import re
import random
import sys
import os
import math
import logging

from subprocess import check_output
from collections import namedtuple, Counter

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbcore.io.GffIO import GffReader, Gff3Record, GffWriter

from pbhla.io.BasH5IO import BasH5Extractor
from pbhla.io.SamIO import SamReader
from pbhla.io.FastaIO import FastaReader, FastaWriter
from pbhla.io.FastqIO import FastqReader, FastqWriter
from pbhla.io.BlasrIO import parse_blasr
from pbhla.io.GffIO import create_annotation, create_var_annotation
from pbhla.stats.AmpliconFinder import AmpliconFinder
from pbhla.stats.SubreadStats import SubreadStats
from pbhla.align.MultiSequenceAligner import MSA_aligner
from pbhla.smrtanalysis.SmrtAnalysisTools import SmrtAnalysisRunner
from pbhla.utils import make_rand_string, getbash, runbash, create_directory
from pbhla.fasta.Utils import write_fasta, fasta_size, extract_sequence, trim_fasta

__version__ = "0.1.0"
 
file_info = namedtuple('file_info', 'fasta_fn, ref_name, locus')
info = namedtuple('info', 'canonical_pos, feature, codon' )

SMRT_ANALYSIS = "/mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh"
NUM_PROC = 4

class HlaPipeline( PBMultiToolRunner ):

    def __init__(self):
        self.parse_args()
        self.initialize_project()
        self.initialize_logging()
        self.validate_settings()

    def parse_args( self ):
        import argparse
        desc = "A pipeline for performing HLA haplotype sequencing."
        super(HlaPipeline, self).__init__(desc)
        parser = argparse.ArgumentParser()
        parser.add_argument("--bash5", 
                help="A list of absolute paths to bas.h5 files containing samples to be HLA typed.")
        parser.add_argument("--proj", 
                help="Indicates that you want to run from an existing HLA typing \n \
                project folder specified here.")
        parser.add_argument("--refs", 
                help="This is a space delimited text file with two columns: (path to fasta) (locus) \n \
                Path to fasta points to a single sequence fasta file. \n \
                Sequences will be used to sort reads into groups, and as a backbone for the phasing of each cluster.")
        parser.add_argument("--choose_ref",
                help="Fofn of multifasta files that contain multiple reference sequences for each locus. \n \
                HLA pipeline will examine the data and choose appropriate references from this multifasta. \n \
                This should be a space delimited text file with two columns: (path to multifasta) (locus)")
        parser.add_argument("--by_locus",
                            help="Separate reads by locus rather than by reference")
        parser.add_argument("--MSA",
                help="This is a fofn which describes which prealigned MSA corresponds to each locus.")
        parser.add_argument("--dilute", dest="dilute_factor", default=str(1.0),
                help="Only use this proportion of the reads from the bas.h5 for phasing. \n \
                It wont affect the resquencing step though, the whole .bas.h5 will still be used there.")
        parser.add_argument("--min_read_score", type=float, default=0.75,
                            help="Only extract subreads more accurate than this from bas.h5 files")
        parser.add_argument("--min_read_length", type=int, default=500,
                            help="Only extract subreads longer than this from bas.h5 files")
        parser.add_argument("--phasr-args", nargs='*', default=[''],
                            help="pass these args to phasr.")
        parser.add_argument("--region_table", 
                            help="Region Table of White-Listed reads to use")
        parser.add_argument("--nproc", type=int, 
                            dest="nproc", default=NUM_PROC,
                            help="Number of processors to use for parallelized computing")
        parser.add_argument("--log", action='store_true',
                            help="Create log file.")
        parser.add_argument("--HGAP_dilute", default=1.0,
                            help="Proportion of reads to pass to HGAP.")
        parser.add_argument("--HGAP", action='store_true',
                            help="Run HGAP on unmapped reads.")
        parser.add_argument("--avoid_phasr", action='store_true',
                            help="Avoid phasr if reference mapping has identified two likely phases.")
        parser.add_argument("--annotate", action='store_true',
                            help="Avoid annotation.")
        self.args = parser.parse_args()

    def initialize_project( self ):
        if self.args.proj == None:
            rand_string=make_rand_string()
            self.args.proj = "%s/results_%s" % (os.getcwd(), rand_string)
            self.welcome_msg="Beginning new HLA calling project. Results will be saved in %s\n" % self.args.proj
            os.mkdir(self.args.proj)
        else:
            if os.path.isdir(str(self.args.proj)):
                self.welcome_msg="Processing existing HLA calling project from %s\n" % self.args.proj
            else:
                try:
                    os.mkdir(self.args.proj)
                    self.welcome_msg="Beginning new HLA calling project. Results will be saved in %s\n" % self.args.proj
                except:
                    raise SystemExit
        # Create the various sub-directories
        self.main_dir = self.args.proj
        self.log_dir = os.path.join( self.main_dir, 'log' )
        create_directory( self.log_dir )
        self.subread_dir = os.path.join( self.main_dir, 'subreads' )
        create_directory( self.subread_dir )
        self.stats_dir = os.path.join( self.main_dir, 'statistics' )
        create_directory( self.stats_dir )
        self.phased_dir = os.path.join( self.main_dir, 'phased' )
        create_directory( self.phased_dir )
        self.reseq_dir = os.path.join( self.main_dir, 'reseq' )
        create_directory( self.reseq_dir )
    
    def initialize_logging( self ):
        time_format = "%Y-%m-%d %I:%M:%S"
        log_format = "%(asctime)s %(levelname)s %(processName)s " + \
                     "%(funcName)s %(lineno)d %(message)s"
        out_format = "%(asctime)s %(funcName)s %(message)s"
        self.logger = logging.getLogger()
        self.logger.setLevel( logging.INFO )
        # Set-up one logger for STDOUT
        h1 = logging.StreamHandler( stream=sys.stdout )
        h1.setLevel( logging.INFO )
        f1 = logging.Formatter( fmt=out_format, datefmt=time_format )
        h1.setFormatter( f1 )
        self.logger.addHandler( h1 )
        # If specified, setup a second logger to log to file
        if self.args.log:
            h2 = logging.FileHandler( self.args.proj+"/hla_pipeline.log" )
            f2 = logging.Formatter( fmt=log_format, datefmt=time_format )
            h2.setFormatter( f2 )
            h2.setLevel( logging.INFO )
            self.logger.addHandler( h2 )
        # Begin by logging both the invocation and welcome message
        self.logger.info( "hla_pipeline invoked:\n\t%s" % " ".join(sys.argv) )
        self.logger.info( self.welcome_msg )

    def validate_settings( self ):
        # Make filepath's absolute if needed
        self.args.bash5 = os.path.abspath( self.args.bash5 )
        # Check dilution factors
        if float(self.args.dilute_factor) <= float(0.0) or float(self.args.dilute_factor) > float(1.0):
            self.logger.info("Dilute factor must be between 0 and 1")
            raise SystemExit
        # parse phasr args
        self.phasr_argstring = ''
        for argument in self.args.phasr_args:
            if ':' in argument:
                param, value = argument.split(":")
                self.phasr_argstring += '--%s %s ' % (param, value)
            elif 'output' in argument or 'cname' in argument:
                pass    
            else:
                self.phasr_argstring += '--%s ' % argument
        # Initialize the SmrtAnalysisTools
        self.smrt_analysis = SmrtAnalysisRunner( SMRT_ANALYSIS, self.log_dir, self.args.nproc )

    def getVersion(self):
        return __version__

    def __call__( self ):
        self.extract_subreads()
        self.create_locus_reference()
        self.align_subreads()
        self.find_amplicons()
        self.summarize_aligned_subreads()
        #self.phase_subreads()   
        #self.resequence()
        #self.annotate()

    def initialize_process(self, proc_name, output_files):
        self.logger.info('Initializing the "{0}" process'.format(proc_name))
        for filename in output_files:
            pass

    def extract_subreads( self ):
        print self.args.bash5
        self.logger.info('Beginning the extraction of subread data')
        # Dump all valid reads from the above files
        self.read_dump_file = self.subread_dir + "/read_dump.fasta"
        if os.path.isfile( self.read_dump_file ):
            self.logger.info('Subread dump-file exists ( %s )' % self.read_dump_file)
            self.logger.info('Skipping subread extraction step')
            return
        self.logger.info('Extracting reads according to the following criteria:')
        self.logger.info('\t\tMinimum Subread Length: %s' % self.args.min_read_length)
        self.logger.info('\t\tMinimum Read Score: %s' % self.args.min_read_score)
        self.logger.info('\t\tSubread Dilution Factor: %s' % self.args.dilute_factor)
        extractor = BasH5Extractor( self.args.bash5, 
                                    self.read_dump_file,
                                    min_length=self.args.min_read_length,
                                    min_score=self.args.min_read_score,
                                    dilution=self.args.dilute_factor)
        self.logger.info('Running the subread extraction process...')
        extractor()
        self.logger.info('Finished the extraction of subread data\n')

    def create_locus_reference( self ):
        # Assign reads to specific loci
        self.logger.info("Creating locus reference file for the data")
        self.locus_dict = {}
        self.reference_sequences = self.args.proj + "/subreads/reference_sequences.fasta"
        self.locus_key = self.args.proj + "/subreads/locus_key.txt"
        # First we check whether
        if os.path.isfile( self.reference_sequences ) and os.path.isfile( self.locus_key ):
            self.logger.info("Found existing locus reference panel ( %s )" % self.locus_key)
            self.logger.info("Reading locus reference panel...")
            with open(self.locus_key, "r") as handle:
                for line in handle:
                    gene, locus = line.strip().split()
                    self.locus_dict[gene] = locus
        # If no locus key and Choose_Ref is None, read the locus from the regular reference
        elif self.args.choose_ref is None:
            self.logger.info("No locus reference panel found, creating one...")
            self.logger.info("Choose-Ref option not specified, using default references")
            ### get list of loci in reference fasta
            with open(self.args.refs, "r") as handle:
                for line in handle:
                    fasta_file, locus = line.strip().split()
                    try:
                        for record in FastaReader( fasta_file ):
                            write_fasta([record], self.reference_sequences, "a")
                            self.locus_dict[record.name] = locus
                            with open(self.locus_key, "a") as of:
                                print>>of, "%s %s" % (record.name, locus)
                            break
                    except:
                        msg = 'Could not read reference fasta "%s"' % fasta_file
                        self.logger.info( msg )
                        raise IOError( msg )
        elif self.args.by_locus:
            raise IOError("Not Yet Implemented")
        # If no locus key and Choose_Ref is set, read the loci from the Choose_Ref file
        else:
            ### if the --choose-ref option is selected we will choose two reference sequences
            ### from the multifasta of sanger references for each locus, based on which two sequences aligned the most reads
            ### we will only choose two references if the two sequences that aligned the most
            ### reads are significantly different, this prevents us from splitting up 
            ### our CLRs if the sample is really homozygous at that locus.
            # TODO: Refactor this block
            with open(self.args.choose_ref, "r") as handle:          
                for line in handle:
                    reffasta, locus = line.strip().split()
                    if os.path.isfile( reffasta ):
                        pass
                    else:
                        self.logger.info("Cound not open multifasta ( %s )" % ( reffasta ) )
                    nCandidates = fasta_size(reffasta)
                    blasr_output = self.args.proj+"/subreads/"+reffasta.split("/")[-1].split(".")[0]+".blasr"
                    if not os.path.isfile(blasr_output):
                        self.logger.info("Aligning all raw reads to the reference multifasta ( %s )." % ( reffasta ) )
                        runbash("blasr %s %s -bestn 1 -nCandidates %s -maxScore -10000 -nproc %s -out %s" % \
                            ( self.read_dump_file, reffasta, nCandidates, self.args.nproc, blasr_output) )
                        if os.path.getsize(blasr_output) <= 1:          
                            self.logger.info("No raw reads aligned suitably to any reference from ( %s )." % (reffasta) )
                            continue
                    else:
                        self.logger.info("Already found blasr output ( %s )." % (blasr_output) )
                        if os.path.getsize(blasr_output) <= 1:
                            self.logger.info("No raw reads aligned suitably to any reference from ( %s )." % (reffasta) )
                            continue
                    ### This unix command will just read how many times a read aligned to each reference
                    ref_align_counts = [ x.split() for x in getbash("awk \'{ print $2 }\' %s | sort | uniq -c | sort -k 1 -g -r" \
                                         % ( blasr_output) ).strip().split("\n")]
                    bestrefs=[]; shift=1; finish_up=False; similar_refs=False
                    best_ref = ref_align_counts.pop(0)
                    sequence = extract_sequence(reffasta, [ best_ref[1] ])
                    write_fasta(sequence, self.reference_sequences, "a")
                    self.locus_dict[sequence[0].name] = locus
                    with open(self.locus_key, "a") as of:
                        print>>of, "%s %s" % (sequence[0].name, locus)
                    ref_fn = self.args.proj+"/subreads/"+reffasta.split("/")[-1]+"_ref"
                    self.logger.info("Best ref ( %s ) aligned ( %s ) reads at score < -10000." % (best_ref[1], best_ref[0]) )
                    while True:
                        try:
                            potential_ref = ref_align_counts.pop(0)
                            self.logger.info("( %s ) aligned ( %s ) reads at score < -10000." % (potential_ref[1], potential_ref[0]) ) 
                        except IndexError:
                            break
                        sequence_pair = extract_sequence(reffasta, [ best_ref[1], potential_ref[1] ] )
                        write_fasta([sequence_pair[0]], ref_fn+"1.fasta", mode="w")
                        write_fasta([sequence_pair[1]], ref_fn+"2.fasta", mode="w")
                        try:
                            pctid = getbash("blasr %s %s -bestn 1" % ( ref_fn+"1.fasta", ref_fn+"2.fasta") ).strip().split()[5]
                            pctid = float(pctid)
                        except:
                            pctid = float(0.0)
                        if float(pctid) > 99.0:
                            self.logger.info("Two refs are very similar. Will pick next best. Pctid ( %s ) > ( 99.0 )" % (pctid) )  
                            continue
                        elif float(pctid) < 99.0:
                            self.logger.info("( %s ) chosen." % (potential_ref[1]) )
                            second_best_ref = potential_ref
                            sequence = extract_sequence(reffasta, [second_best_ref[1] ] )
                            write_fasta(sequence, self.reference_sequences, "a")
                            self.locus_dict[sequence[0].name] = locus
                            with open(self.locus_key, "a") as of:
                                print>>of, "%s %s" % (sequence[0].name, locus)
                            break
        self.logger.info("Finished creating locus reference panel\n")

    def align_subreads( self ):
        self.logger.info("Beginning alignment of data to selected references")
        self.sam_file = self.args.proj+"/subreads/read_dump.sam"   
        if os.path.isfile( self.sam_file ):
            self.logger.info("Found existing SAM file ( %s )" % self.sam_file)
            self.logger.info("Skipping alignment step...\n")
            return
        ref_count = fasta_size( self.reference_sequences )
        self.logger.info( "Aligning subreads to ( %s ) reference sequences" % ref_count)
        # TODO: Move Blasr calls to their own module
        blasr_cline = "blasr %s %s -bestn 1 -nCandidates %s -noSplitSubreads -sam -nproc %s -out %s" % (self.read_dump_file, 
                                                                                                        self.reference_sequences, 
                                                                                                        ref_count, 
                                                                                                        self.args.nproc, 
                                                                                                        self.sam_file)
        runbash( blasr_cline )
        try:
            assert os.path.isfile( self.sam_file )
        except:
            msg = 'Blasr returned no output for locus assignment'
            self.logger.info( msg )
            raise IOError( msg )
        self.logger.info("Finished aligning data to selected references\n")

    def find_amplicons(self):
        self.logger.info("Finding amplicon locations within each reference")
        self.amplicon_summary = self.args.proj + "/subreads/amplicon_summary.csv"
        if os.path.isfile( self.amplicon_summary ):
            self.logger.info("Already found Amplicon Summary ( {0} )".format(self.amplicon_summary))
            self.logger.info("Skipping amplicon finding step...\n")
            return
        finder = AmpliconFinder(self.sam_file, self.amplicon_summary)
        finder()
        self.logger.info("Finished finding amplicon locations\n")

    def summarize_aligned_subreads(self):
        self.logger.info("Assigning subreads to their associated locus")
        self.subread_files = self.args.proj + "/subreads/subread_files.txt"
        self.unmapped_reads = self.args.proj + "/subreads/unmapped_reads.fasta"
        self.locus_stats = os.path.join(self.stats_dir, "locus_statistics.csv")
        self.reference_stats = os.path.join(self.stats_dir, "reference_statistics.csv")
        self.amplicon_stats = os.path.join(self.stats_dir, "amplicon_statistics.csv")
        ### will go through sam file, sort out the raw reads, and tabulate statistics while we do it
        if os.path.isfile( self.subread_files ):
            self.logger.info("Already found subread files ( %s )" % self.subread_files )
            self.logger.info("Skipping locus assignment step...\n")
            return
        read_cluster_map = {}
        cluster_files = {}
        stats = SubreadStats( self.reference_sequences, self.locus_dict, self.amplicon_summary )
        # First
        for alignment in SamReader( self.sam_file ):
            read_cluster_map[alignment.qname] = alignment.rname
            locus = self.locus_dict[alignment.rname]
            stats.add_aligned_read( locus, alignment )
        # Next
        for record in FastaReader( self.read_dump_file ):
            try:
                name = read_cluster_map[record.name]
                locus = self.locus_dict[read_cluster_map[record.name]]
                cluster_file = self.args.proj + "/subreads/" + name + "_mapped_reads.fasta"
                write_fasta([record], cluster_file, "a")
                cluster_files[name] = [ cluster_file, name, locus ]
            except KeyError:
                write_fasta([record], self.unmapped_reads, "a")
        ### write out a text file that can act a legend for which fasta file of clustered reads
        ### corresponds to which template from the reference file and which locus
        with open( self.subread_files , "w") as of:
            for item in cluster_files.itervalues():
                print >>of, "%s %s %s" % ( item[0], item[1], item[2]  )
            print >>of, "%s %s %s" % ( self.unmapped_reads, "unmapped", "unmapped") 
        ### finally write out the subread statistics
        stats.write( self.locus_stats, 'locus' )
        stats.write( self.reference_stats, 'reference' )
        stats.write( self.amplicon_stats, 'amplicon' )
        self.logger.info("Subread files and summary created ( %s )\n" % self.subread_files )

    def phase_subreads( self ):
        self.logger.info("Starting the sub-read phasing process")
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
                self.logger.info( msg )
                raise IOError( msg )

            if os.path.isfile(phasr_output):
                for record in FastaReader(phasr_output):
                    seq_to_locus[record.name] = locus
                self.logger.info("Found phasr output for ( %s ) in ( %s ). Skipping." % (fasta_file, phasr_output))
                continue
            
            ### if this locus is not to be phased, we can just set max recursion level to 0, and phasr will build regular consensus
            if self.args.avoid_phasr and not to_be_phased[locus]:
                self.logger.info('Running Gcon instead of Phasr for "%s"' % fasta_file)
                argstring = re.sub('--max_recursion_level [0-9]',
                                   '--max_recursion_level 0',
                                   self.phasr_argstring)
                argstring += ' --max_recursion_level 0'
            else:
                self.logger.info('Running  Phasr on "%s"' % fasta_file)
                argstring = self.phasr_argstring
            phasr_cmd = "phasr %s --ref %s\
                --output %s \
                --log \
                %s" % ( fasta_file, backbone_fasta, phasr_output, self.phasr_argstring)
            ### run phasr
            self.logger.info("phasr invoked: %s" % phasr_cmd)
            check_output(". /home/UNIXHOME/jquinn/HGAP_env/bin/activate; %s" % (phasr_cmd),  
                         executable='/bin/bash', shell=True)
            if os.path.isfile( phasr_output ):
                output_count = fasta_size(phasr_output)
                self.logger.info('Phasr output %s seqs to "%s"' % ( output_count, phasr_output) )
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
        self.logger.info("Phasr created %s sequence(s) total" % phasr_output_count )
        self.logger.info("Phasing complete.\n")

    def resequence(self, step=1):
        ### TODO: maybe adjust compare seq options so that euqal mapping are not assigned randomly...
        ### we resequence once, elminate redundant clusters, then resequence the remaining clusters
        ### to avoid splitting our CLRs over clusters that respresent the same template 
        self.logger.info("Resequencing step %s begun." % ( str(step) ))
        reseq_dir = os.path.join( self.args.proj + '/reseq' )
        create_directory( reseq_dir )
        
        reseq_references = self.args.proj+"/reseq/reseq_references_%s.fasta" % step
        reseq_output = self.args.proj+"/reseq/reseq_output_%s" % step
        reseq_fasta = reseq_output + ".fasta"
        reseq_fastq = reseq_output + ".fastq"
        reseq_variant_gff = reseq_output + "_variants.gff"
        reseq_coverage_gff = reseq_output + "_coverage.gff"
        reseq_coverage = reseq_output + "_coverage.txt"

        #if os.path.isfile( reseq_output+".fasta" ):
        #    self.logger.info("Already found reseq output at ( %s )\n" % ( reseq_output+".fasta") )
        #    return

        if not os.path.isfile( reseq_references ):
            if step == 1:
                ### combine the denovo preassembled trimmed reads with the phasr output
                try:
                    runbash(" cat %s %s > %s" % (self.args.proj+"/unmapped/pre_assembled_reads_trimmed.fasta", self.hap_con_fasta, reseq_references) ) 
                except:
                    runbash(" cat %s > %s" % (self.hap_con_fasta, reseq_references) )

        ### create a .cmp.h5
        self.reference_alignment = os.path.join( reseq_dir, "reference.cmp.h5")
        if os.path.isfile( self.reference_alignment ):
            self.logger.info("Found existing Reference Alignment ( %s )" % self.reference_alignment)
            self.logger.info("Skipping Reference Alignment step...\n")
        else:
            self.smrt_analysis.compare_sequences( self.args.bash5, 
                                                  reseq_references,
                                                  self.reference_alignment )
            self.smrt_analysis.load_pulses( self.args.bash5, self.reference_alignment )
            self.smrt_analysis.sort_cmpH5( self.reference_alignment )
        
        ### run quiver 
        if os.path.isfile( reseq_fasta ):
            self.logger.info("Found existing Quiver Results ( %s )" % reseq_fasta)
            self.logger.info("Skipping Quiver Resequencing step...\n")
        else:
            self.smrt_analysis.quiver( self.reference_alignment,
                                       reseq_references, 
                                       reseq_fasta, 
                                       reseq_fastq,
                                       reseq_variant_gff)
        
        ## Summarize the coverage
        self.consensus_alignment = os.path.join( reseq_dir, "consensus.cmp.h5" )
        if os.path.isfile( self.consensus_alignment ):
            self.logger.info("Found existing Consensus Alignment ( %s )" % self.consensus_alignment)
            self.logger.info("Skipping Consensus Alignment step...\n")
        else:
            self.smrt_analysis.compare_sequences( self.args.bash5,
                                                  reseq_fasta,
                                                  self.consensus_alignment)
            self.smrt_analysis.sort_cmpH5( self.consensus_alignment )

        ref_dir = "reseq_references_" + str(step)
        ref_dir_path = os.path.join(self.args.proj, 'reseq', ref_dir)
        if os.path.isdir( ref_dir_path ):
            self.logger.info("Found existing Reference Directory ( %s )" % ref_dir_path)
            self.logger.info("Skipping Reference Creation step...\n")
        else:
            self.smrt_analysis.reference_uploader( reseq_fasta, ref_dir, reseq_dir ) 
       
        if os.path.isfile( reseq_coverage_gff ):
            self.logger.info("Found existing Coverage Analysis file ( %s )" % reseq_coverage_gff)
            self.logger.info("Skipping Coverage Analysis step...\n")
        else:
            self.smrt_analysis.summarize_coverage( self.consensus_alignment,
                                                   ref_dir_path, 
                                                   reseq_coverage_gff ) 

        if os.path.isfile( reseq_coverage ):
            self.logger.info("Found existing Coverage Summary file ( %s )" % reseq_coverage)
            self.logger.info("Skipping Coverage Summary step...\n")
        else:
            runbash("sed \'/^#/d\' %s | awk \'function getcov (entry) { split(entry,a,\",\"); out=substr(a[1],6,length(a[1])); return out } \
                BEGIN{ n=0; t=0; print \"Seq_Name\",\"Avg_Cov\" }\
                { if( NR == 1) { n+=1; t+=getcov($9); last=$1 } else { if( $1 == last ) { n+=1; t+=getcov($9) } else { print last,(t/n) ; n=0; t=0; n+=1; t+=getcov($9); last=$1 } } }\
                END{ print last,(t/n) }\' > %s" % ( reseq_coverage_gff, 
                                                    reseq_coverage ))

        reseq_quality = self.args.proj + "/reseq/quality_summary_%s.txt" % step
        if os.path.isfile( reseq_quality ):
            self.logger.info("Found existing Quality Summary file ( %s )" % reseq_quality)
            self.logger.info("Skipping Quality Summary step...\n")
        else:
            with open( reseq_quality, "w" ) as handle:
                for record in FastqReader( reseq_fastq ):
                    average_quality = sum(record.quality)/float(len(record.quality))
                    handle.write("%s %s\n" % (record.name, average_quality))

        ### separate out resequenced HLA sequences from the nonspecific ones, by sequence name, store in two different files
        ### TODO: read fastq output, parse consnesus QV, store somewhere the annotator can get at it
        target_seq_names=[]
        for record in FastaReader(self.hap_con_fasta):
            target_seq_names.append( record.name )
        with FastaWriter(self.args.proj+"/reseq/resequenced_hap_con.fasta") as on_target: 
            with FastaWriter(self.args.proj+"/reseq/resequenced_nonspecific.fasta") as off_target:
                for record in FastaReader(reseq_output+".fasta"):
                    if record.name.split("|")[0] in target_seq_names:
                        on_target.writeRecord( record )
                    else:
                        off_target.writeRecord( record )
        with FastqWriter(self.args.proj+"/reseq/resequenced_hap_con.fastq") as on_target:
            with FastqWriter(self.args.proj+"/reseq/resequenced_nonspecific.fastq") as off_target:
                for record in FastqReader(reseq_output+".fastq"):
                    if record.name.split("|")[0] in target_seq_names:
                        on_target.writeRecord( record )
                    else:   
                        off_target.writeRecord( record )

    def annotate(self):     
        if not self.args.annotate:
            return
        try:
            os.mkdir(self.args.proj+"/annotate")
        except:
            pass
        self.logger.info("Annotation begun")
        self.locus_dict={}       
        nseqs=0
        with open(self.args.proj+"/phased/phasr_output_seqs.txt", "r") as f:
            for line in f:
                self.locus_dict[line.strip().split()[0]] = line.strip().split()[1]
                nseqs+=1
        self.logger.info("( %s ) sequences to annotate." % (nseqs) )
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
            self.logger.info("Processing ( %s ) from locus ( %s )." % (name, locus) )

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
                    self.logger.info("MSA failed for ( %s )" % (name) )
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
                    self.logger.info("MSA failed for ( %s )" % (name+"_cDNA") )
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
                    self.logger.info("MSA failed for ( %s )" % (locus) )
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
        self.logger.info("There are ( %s ) unmapped reads, ( %s ) percent will be passed to HGAP." % ( num_unmapped_reads, self.args.HGAP_dilute) )
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
        self.logger.info("HGAP envoked")
        output = check_output( "source /home/UNIXHOME/jchin/Share/HGA_env/bin/activate; cd %s/unmapped; %s %s" % (self.args.proj, HGAP_executable, self.args.proj+"/unmapped/HGAP.cfg"), executable='/bin/bash', shell=True )       
        self.logger.info("HGAP exited with the following output:\n%s" % output )
        runbash("cat %s > %s" % (self.args.proj+"/unmapped/dist_map/pre_assembled_reads_*.fa", self.args.proj+"/unmapped/pre_assembled_reads.fasta") )
        self.logger.info("HGAP created ( %s ) PLRs." % ( fasta_size(self.args.proj+"/unmapped/pre_assembled_reads.fasta") ) )
        trim_fasta(self.args.proj+"/unmapped/pre_assembled_reads.fasta", 
            self.args.proj+"/unmapped/pre_assembled_reads_trimmed.fasta", 
            mode="redundant", 
            pctid_threshold=float(90), 
            coverage=None, 
            aln_portion=float(0.50) )   
        self.logger.info("After trimming there are ( %s ) unique PLRs" % ( fasta_size(self.args.proj+"/unmapped/pre_assembled_reads_trimmed.fasta") ) )
"""
if __name__ == '__main__':
    pipeline = HlaPipeline()
    pipeline()
