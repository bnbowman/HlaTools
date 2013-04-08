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

from pbtools.pbhla.io.SamIO import SamReader
from pbtools.pbhla.io.FastaIO import FastaReader, FastaWriter
from pbtools.pbhla.io.FastqIO import FastqReader, FastqWriter
from pbtools.pbhla.io.BlasrIO import parse_blasr
from pbtools.pbhla.io.GffIO import create_annotation, create_var_annotation
from pbtools.pbhla.stats.SubreadStats import SubreadStats
from pbtools.pbhla.align.MultiSequenceAligner import MSA_aligner
from pbtools.pbhla.utils import make_rand_string, getbash, runbash
from pbtools.pbhla.fasta.Utils import write_fasta, fasta_size, extract_sequence, trim_fasta

__version__ = "0.1.0"
 
class HLApipeline(PBMultiToolRunner):

	def __init__(self):
                import argparse
		desc = "A pipeline for performing HLA haplotype sequencing."
		super(HLApipeline, self).__init__(desc)
		args = argparse.ArgumentParser()
		args.add_argument("--bash5", dest="bash5",
			    help="A list of absolute paths to bas.h5 files containing samples to be HLA typed.")
		args.add_argument("--proj", dest="proj",
			    help="Indicates that you want to run from an existing HLA typing \n \
			    project folder specified here.")
		args.add_argument("--refs", dest="refs",
			    help="This is a space delimited text file with two columns: (path to fasta) (locus) \n \
				Path to fasta points to a single sequence fasta file. \n \
				Sequences will be used to sort reads into groups, and as a backbone for the phasing of each cluster.")
		args.add_argument("--choose-ref", dest="choose_ref", default=None,
			    help="Fofn of multifasta files that contain multiple reference sequences for each locus. \n \
				HLA pipeline will examine the data and choose appropriate references from this multifasta. \n \
                                This should be a space delimited text file with two columns: (path to multifasta) (locus)")
		args.add_argument("--MSA", dest="MSA",
			    help="This is a fofn which describes which prealigned MSA corresponds to each locus.")
		args.add_argument("--dilute", dest="dilute_factor", default=str(1.0),
			    help="Only use this proportion of the reads from the bas.h5 for phasing. \n \
				It wont affect the resquencing step though, the whole .bas.h5 will still be used there.")
		args.add_argument("--min-read-length", dest="min_read_length", default=str(500),
			    help="Only extract reads longer than this from the bas.h5.")	
		args.add_argument("--phasr-args", dest="phasr_args", nargs='*', default=[''],
			    help="pass these args to phasr.")
		args.add_argument("--nproc", dest="nproc", default=int(1),
                            help="Number of processors to use for parallelized computing")
		args.add_argument("--log", dest="log", action='store_true', default=False,
			    help="Create log file.")
		args.add_argument("--HGAP_dilute", dest="HGAP_dilute", default=int(1),
                            help="Proportion of reads to pass to HGAP.")
		args.add_argument("--HGAP", dest="HGAP", default=False,
                            help="Run HGAP on unmapped reads.")
		args.add_argument("--avoid_phasr", dest="avoid_phasr", default=False,
                            help="Avoid phasr if reference mapping has identified two likely phases.")
		args.add_argument("--annotate", dest="annotate", action = 'store_true', default=False,
                            help="Avoid annotation.")
		self.args = args.parse_args()


		if self.args.proj == None:
			rand_string=make_rand_string()
			self.args.proj = "%s/results_%s" % (os.getcwd(), rand_string)
			welcome_msg="Beginning new HLA calling project. Results will be saved in %s" % self.args.proj
			os.mkdir(self.args.proj)
		else:
			if os.path.isdir(str(self.args.proj)):
				welcome_msg="Processing existing HLA calling project from %s" % self.args.proj
			else:
				try:
					os.mkdir(self.args.proj)
					welcome_msg="Beginning new HLA calling project. Results will be saved in %s" % self.args.proj
				except:
					sys.exit(0)
	
	        self.logger = logging.getLogger()
                self.logger.setLevel(logging.INFO)
                if self.args.log:
                    h1 = logging.FileHandler(self.args.proj+"/hla_pipeline.log")
                elif not self.args.log:
                    h1 = logging.FileHandler("/dev/null")
                f = logging.Formatter("%(processName)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
                h1.setFormatter(f)
                h1.setLevel(logging.INFO)
                self.logger.addHandler(h1)
		h2 = logging.StreamHandler(stream=sys.stdout)
		h2.setFormatter(f); h2.setLevel(logging.INFO); self.logger.addHandler(h2)
                self.logger.info("hla_pipeline envoked: %s" % " ".join(sys.argv) )
		self.logger.info(welcome_msg)

		if float(self.args.dilute_factor) <= float(0.0) or float(self.args.dilute_factor) > float(1.0):
			self.logger.info("Dilute factor must be between 0 and 1")
			sys.exit(0)

	def getVersion(self):
		return __version__

	def subread_extraction(self):
	    try:
		os.mkdir(self.args.proj+"/subreads")
	    except:
		pass

	    if not os.path.isfile(str(self.args.bash5)):
		self.logger.info("Could not open ( %s )." % self.args.bash5 ) 	
	    else:
		### parse bash5 fofn
		with open(self.args.bash5) as f:
		    bash5_files=[]
		    for line in f.readlines():
			    if os.path.isfile(line.strip()): 
				bash5_files.append(line.strip())
			    else:
				self.logger.info("Could not open ( %s )." % line.strip() )
		self.logger.info("Subread extraction begun on the following files:\n%s" % '\n'.join(bash5_files))

	    read_dump_fn=self.args.proj+"/subreads/read_dump.fasta"
	    if os.path.isfile(read_dump_fn):
		self.logger.info("Already found subread dump for file ( %s )" % f)
	    else:
		for f in bash5_files:
		    ### dump all reads > min_readscore into a fasta
		    ### but check the project directory to see if this has been done already
		    ### and dont repeat if the read dump is already there
		    #self.logger.info( "Dumping reads from ( %s ) with score > ( %s )" % ( f, self.args.min_readscore ) )
		    self.logger.info( "Dumping reads from ( %s ) with length > ( %s )" % ( f, self.args.min_read_length ) )
		    runbash("dumpReads.py "+f+" "+self.args.dilute_factor+" "+self.args.min_read_length+" >> "+read_dump_fn )

	    self.logger.info("Assuming all bash5 files belong to the same sample, reads will be assigned to each locus based on best alignment..")
	    locus_dict={}
	    if self.args.choose_ref == None and not os.path.isfile(self.args.proj+"/subreads/reference_sequences.fasta"):
		### get list of loci in reference fasta
		with open(self.args.refs, "r") as f:
		    for line in f:
			fasta_fn = line.split()[0]
			locus = line.split()[1]
			try:
			    f2 = FastaReader(fasta_fn)
			    for r in f2:
				write_fasta([r], self.args.proj+"/subreads/reference_sequences.fasta", "a")
				locus_dict[r.name]=locus
				with open(self.args.proj+"/subreads/locus_key.txt", "a") as of:
				    print>>of, "%s %s" % (r.name, locus)
				break
			except:
			    self.logger.info("Could not read from reference fasta ( %s )" % fasta_fn )
	    
	    ### bin out reads based on which reference they align to
	    if self.args.choose_ref != None and not os.path.isfile(self.args.proj+"/subreads/reference_sequences.fasta"):
		with open(self.args.choose_ref, "r") as f:		
		    ### if the --choose-ref option is selected we will choose two reference sequences
		    ### from the multifasta of sanger references for each locus, based on which two sequences aligned the most reads
		    ### we will only choose two references if the two sequences that aligned the most
		    ### reads are significantly different, this prevents us from splitting up 
		    ### our CLRs if the sample is really homozygous at that locus.
		    for line in f:
			reffasta = line.strip().split()[0]
			locus = line.strip().split()[1]
			if os.path.isfile(reffasta):
			    pass
			else:
			    self.logger.info("Cound not open multifasta ( %s )" % ( reffasta ) )
			    continue
			
			nCandidates = fasta_size(reffasta)
			blasr_output = self.args.proj+"/subreads/"+reffasta.split("/")[-1].split(".")[0]+".blasr"
			if not os.path.isfile(blasr_output):
			    self.logger.info("Aligning all raw reads to the reference multifasta ( %s )." % ( reffasta ) )
			    runbash("blasr %s %s -bestn 1 -nCandidates %s -maxScore -10000 -nproc %s -out %s" % \
				( read_dump_fn, reffasta, nCandidates, self.args.nproc, blasr_output) )
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
			write_fasta(sequence, self.args.proj+"/subreads/reference_sequences.fasta", "a")
			locus_dict[sequence[0].name]=locus
			with open(self.args.proj+"/subreads/locus_key.txt", "a") as of:
			    print>>of, "%s %s" % (sequence[0].name, locus)
			ref_fn = self.args.proj+"/subreads/"+reffasta.split("/")[-1]+"_ref"
			self.logger.info("Best ref ( %s ) aligned ( %s ) reads at score < -10000." % (best_ref[1], best_ref[0]) )
			while 1:
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
				write_fasta(sequence, self.args.proj+"/subreads/reference_sequences.fasta", "a")
				locus_dict[sequence[0].name] = locus
				with open(self.args.proj+"/subreads/locus_key.txt", "a") as of:
				    print>>of, "%s %s" % (sequence[0].name, locus)
				break
			

	    if os.path.isfile(self.args.proj+"/subreads/reference_sequences.fasta") and os.path.isfile(self.args.proj+"/subreads/locus_key.txt"):
		with open(self.args.proj+"/subreads/locus_key.txt", "r") as f:
		    for line in f:
			locus_dict[line.strip().split()[0]] = line.strip().split()[1]
		self.logger.info("Already found reference panel.")
		    
	    sam_fn = self.args.proj+"/subreads/read_dump.sam"	
	    if not os.path.isfile(sam_fn):
		f = FastaReader(self.args.proj+"/subreads/reference_sequences.fasta")
		self.logger.info( "Aligning all raw reads to ( %s ) reference sequences" % ( fasta_size(self.args.proj+"/subreads/reference_sequences.fasta") ) )
		runbash("blasr %s %s -bestn 1 -nCandidates %s -noSplitSubreads -sam -nproc %s -out %s" % \
		    (read_dump_fn, self.args.proj+"/subreads/reference_sequences.fasta", fasta_size(self.args.proj+"/subreads/reference_sequences.fasta"), self.args.nproc, sam_fn) )
		if not os.path.isfile(sam_fn):
		    self.logger.info("Blasr returned no output for locus assignment. Exiting subread extraction step.")
		    return 0
	    else:
		self.logger.info("Already found sam file at ( %s ). Will skip locus assignment." % sam_fn)

	    ### will go through sam file, sort out the raw reads, and tabulate statistics while we do it
	    if not os.path.isfile(self.args.proj+"/subreads/subread_files.txt"):
		read_cluster_map={}; subread_files={}
		stats = SubreadStats(self.args.proj+"/subreads/reference_sequences.fasta", locus_dict)
		f = SamReader(sam_fn)
		for alignment in f:
		    stats.loci[locus_dict[alignment.rname]]['aligned_reads']+=1
		    if 'fp' in alignment.qname: stats.loci[locus_dict[alignment.rname]]['aligned_fp_reads']+=1
		    stats.loci[locus_dict[alignment.rname]]['aligned_bp']+=alignment.aln_length
		    if 'fp' in alignment.qname: stats.loci[locus_dict[alignment.rname]]['aligned_bp']+=alignment.aln_length
		    read_cluster_map[alignment.qname]=alignment.rname
			
		f = FastaReader(read_dump_fn)
		for r in f:
		    stats.total_reads+=1
		    if 'fp' in r.name: stats.total_fp_reads+=1
		    try:
			name=read_cluster_map[r.name]
			locus=locus_dict[read_cluster_map[r.name]]
			stats.loci[locus]['raw_bp']+=len(r.sequence)
			if 'fp' in r.name: stats.loci[locus]['raw_fp_bp']+=len(r.sequence)
			write_fasta([r], self.args.proj+"/subreads/"+name+"_mapped_reads.fasta", "a")
			subread_files[name] = [ self.args.proj+"/subreads/"+name+"_mapped_reads.fasta", name, locus]
		    except KeyError:
			write_fasta([r], self.args.proj+"/subreads/unmapped_reads.fasta", "a")
		### write out a text file that can act a legend for which fasta file of clustered reads
		### corresponds to which template from the reference file and which locus
		with open(self.args.proj+"/subreads/subread_files.txt", "w") as of:
		    for item in subread_files.itervalues():
			print >>of, "%s %s %s" % ( item[0], item[1], item[2]  )
		    print >>of, "%s %s %s" % ( self.args.proj+"/subreads/unmapped_reads.fasta", "unmapped", "unmapped") 
		    
		### finally write out the subread statistics
		try:
		    os.mkdir(self.args.proj+"/statistics")
		except:
		    pass
		stats.write(self.args.proj+"/statistics/")
	    else:
		self.logger.info("Already found subread files at ( %s )" % ( self.args.proj+"/subreads/subread_files.txt" ) )
	    return 0
		
	def phasing(self):
	    if not os.path.isfile(self.args.proj+"/subreads/subread_files.txt"):
		self.logger.info("Phasing skipped.")
		return 0    
            self.logger.info("Phasing begun.")
            try:
                    os.mkdir(self.args.proj+"/phased")
            except:
                    pass
            with open(self.args.proj+"/subreads/subread_files.txt") as f:
                    subread_files=[]
                    file_info = namedtuple('file_info', 'fasta_fn, ref_name, locus')
                    for line in f:
                        subread_files.append(file_info._make(line.strip().split()))
            phasr_output_files=[]

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
            
            for item in subread_files:
                ref_name = item.ref_name
                locus = item.locus
                fasta_fn = item.fasta_fn
		### if these are the unmapped reads we are going to run them through HGAp
		if locus == "unmapped":
		    if not self.args.HGAP:
			continue
		    try:
			os.mkdir(self.args.proj+"/unmapped")	
		    except:
			pass
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
"""	% ( self.args.proj+"/unmapped/input.fofn" ) 	
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
		    continue
				
			
                backbone_fasta=self.args.proj+"/phased/backbone.fa"
                phasr_output=self.args.proj+"/phased/"+ref_name+"_haplotype_consensus.fa"

                ### extract the reference sequence that was used to "cluster" this group of reads
                sequence = extract_sequence(self.args.proj+"/subreads/reference_sequences.fasta", [ ref_name ] )
                write_fasta(sequence, backbone_fasta, "w")
		
                	
                
                if not os.path.isfile(fasta_fn):
                    self.logger.info("Could not open ( %s ). Skipping." % (fasta_fn) )
                    continue
                if os.path.isfile(phasr_output):
		    f = FastaReader(phasr_output)
                    for r in f:
                        seq_to_locus[r.name] = locus
                    self.logger.info("Already found phasr output for ( %s ) in ( %s ). Skipping." % (fasta_fn, phasr_output))
                    continue
                
                ### parse phasr args
                argstring=''
                for argument in self.args.phasr_args:
                        if ':' in argument:
                            argument = argument.split(":")
                            argstring+='--'+argument[0]+' '+argument[1]+' '
                        elif 'output' in argument or 'cname' in argument:
                            pass	
                        else:
                            argstring+='--'+argument+' '
                ### if this locus is not to be phased, we can just set max recursion level to 0, and phasr will build regular consensus
		if self.args.avoid_phasr:
		    if not to_be_phased[locus]:
			argstring = re.sub('--max_recursion_level [0-9]','--max_recursion_level 0',argstring)
			argstring+=' --max_recursion_level 0'
			self.logger.info("Phasr will not be run on ( %s ). Will run graph consensus." % fasta_fn )
                phasr_cmd = "phasr %s --ref %s\
                    --output %s \
                    --log \
                    %s" % ( fasta_fn, backbone_fasta, phasr_output, argstring)
                ### run phasr
                self.logger.info("phasr envoked: %s" % phasr_cmd)
                check_output(". /home/UNIXHOME/jquinn/HGAP_env/bin/activate; %s" % (phasr_cmd),  executable='/bin/bash', shell=True)
		if os.path.isfile(phasr_output):
		    self.logger.info("Phasr output ( %s ) seqs to ( %s )" % ( fasta_size(phasr_output), phasr_output) )

		    ### append results to multifasta
		    runbash("cat %s >> %s" % (phasr_output, self.args.proj+"/phased/haplotype_consensus_sequences.fasta") )
		    f = FastaReader(phasr_output)
		    for r in f:
			seq_to_locus[r.name] = locus
			    

		    ### create record of files created
		    with open(self.args.proj+"/phased/phasr_output_files.txt", "a") as of:
			print >>of, "%s %s %s" % (phasr_output, ref_name, locus)
	    with open(self.args.proj+"/phased/phasr_output_seqs.txt", "w") as of:
		for item in seq_to_locus.iteritems():
		    print >>of, "%s %s" % (item[0], item[1])
	    self.logger.info("Phasr created ( %s ) sequence(s) total" % ( fasta_size(self.args.proj+"/phased/haplotype_consensus_sequences.fasta") ) )
		
            self.logger.info("Phasing complete.")
            return 0
	    

	def resequence(self, step=1):
	    ### TODO: maybe adjust compare seq options so that euqal mapping are not assigned randomly...
	    ### we resequence once, elminate redundant clusters, then resequence the remaining clusters
	    ### to avoid splitting our CLRs over clusters that respresent the same template 
	    self.logger.info("Resequencing step %s begun." % ( str(step) ))
	    try:
		os.mkdir(self.args.proj+"/reseq")
            except:
		pass
			
	    reseq_references=self.args.proj+"/reseq/reseq_references_%s.fasta" % (step)
	    reseq_output=self.args.proj+"/reseq/reseq_output_%s" % (step)
	    if os.path.isfile(reseq_output+".fasta"): 
		self.logger.info("Already found reseq output at ( %s )." % ( reseq_output+".fasta") )
	    if not os.path.isfile(reseq_references):
		if step == 1:
		    ### combine the denovo preassembled trimmed reads with the phasr output
		    try:
			runbash(" cat %s %s > %s" % (self.args.proj+"/unmapped/pre_assembled_reads_trimmed.fasta", self.args.proj+"/phased/haplotype_consensus_sequences.fasta", reseq_references) ) 
		    except:
			runbash(" cat %s > %s" % (self.args.proj+"/phased/haplotype_consensus_sequences.fasta", reseq_references) )

	    ### create a .cmp.h5
	    if not os.path.isfile(self.args.proj+"/reseq/reads.cmp.h5"):
		self.logger.info("compareSequences.py begun.")
		output = check_output("source /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
		compareSequences.py \
			--minAccuracy=0.75 \
			--minLength=50 \
			--placeRepeatsRandomly \
			--minLength=500 \
			--useQuality \
			--h5pbi \
			--info \
			--nproc=4 \
			-x -bestn 1 -x -nCandidates 30 \
			--refSeqName=HLA \
			--algorithm=blasr \
			--noiseData=\"-77.27,0.08654,0.00121\" \
			--h5fn=%s  \
			%s \
			%s" % ( self.args.proj+"/reseq/reads.cmp.h5", 
				self.args.bash5, 
				reseq_references ), executable='/bin/bash', shell=True)	
		self.logger.info("Compare sequences executed with the following messages:\n%s" % output)

		### load pulses to .cmp.h5
		self.logger.info("Load pulses begun.")
		output = check_output("source /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
		loadPulses \
		%s \
		%s \
		-metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag" % \
		( self.args.bash5, self.args.proj+"/reseq/reads.cmp.h5" ), executable='/bin/bash', shell=True )
		self.logger.info("Load pulses executed with the following messages:\n%s" % output)

	    if not os.path.isfile(reseq_output+".fastq"):
		### run quiver 
		output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
			cmph5tools.py sort %s; \
			variantCaller.py \
				--noEvidenceConsensusCall reference \
				--algorithm=quiver \
				-j4  \
				-r %s \
				-o %s \
				-o %s \
				-o %s \
				%s " % \
			( self.args.proj+"/reseq/reads.cmp.h5",
			reseq_references, 
			self.args.proj+"/reseq/variants.gff", 
			reseq_output+".fasta", 
			reseq_output+".fastq",
			self.args.proj+"/reseq/reads.cmp.h5" ), executable='/bin/bash', shell=True)
		self.logger.info("Variant caller executed with the following messages:\n%s" % output)

	    if not os.path.isfile(self.args.proj+"/reseq/coverage_summary_%s.txt" % ( step) ):
		### make new cmph5 after this previous seq cleaning step
		self.logger.info("compareSequences.py begun.")
		output = check_output("source /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
		compareSequences.py \
			--minAccuracy=0.75 \
			--placeRepeatsRandomly \
			--minLength=500 \
			--useQuality \
			--h5pbi \
			--info \
			--nproc=4 \
			-x -bestn 1 -x -nCandidates 30 \
			--refSeqName=HLA \
			--algorithm=blasr \
			--noiseData=\"-77.27,0.08654,0.00121\" \
			--h5fn=%s  \
			%s \
			%s" % ( self.args.proj+"/reseq/reads.cmp.h5",
				self.args.bash5,
				reseq_output+".fasta" ), executable='/bin/bash', shell=True)
		self.logger.info("Compare sequences executed with the following messages:\n%s" % output)

		### get coverage info
		self.logger.info("Coverage summary begun.")
		output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
			referenceUploader -c -f %s -n %s -p %s 2>&1" % (reseq_output+".fasta", "reseq_references_"+str(step), self.args.proj+"/reseq/"), 
			executable='/bin/bash', shell=True)
		ref_dir = self.args.proj+"/reseq/reseq_references_"+str(step)
		output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
					cmph5tools.py sort %s; \
		summarizeCoverage.py %s --reference %s --regionSize 1 > %s" % (self.args.proj+"/reseq/reads.cmp.h5", self.args.proj+"/reseq/reads.cmp.h5", ref_dir, self.args.proj+"/reseq/coverage_%s.gff" % ( step) ),
			executable='/bin/bash', shell=True)

		### create a .txt with two space delimited columns, column 1 = seq name, column 2 = avg coverage
		### TODO: maybe make this not be an awk oneliner for readability
		runbash("sed \'/^#/d\' %s | awk \'function getcov (entry) { split(entry,a,\",\"); out=substr(a[1],6,length(a[1])); return out } \
			BEGIN{ n=0; t=0; print \"Seq_Name\",\"Avg_Cov\" }\
			{ if( NR == 1) { n+=1; t+=getcov($9); last=$1 } else { if( $1 == last ) { n+=1; t+=getcov($9) } else { print last,(t/n) ; n=0; t=0; n+=1; t+=getcov($9); last=$1 } } }\
			END{ print last,(t/n) }\' > %s" % ( self.args.proj+"/reseq/coverage_%s.gff" % (step), self.args.proj+"/reseq/coverage_summary_%s.txt" % (step) ) )

	    if not os.path.isfile(self.args.proj+"/reseq/quality_summary_%s.txt" % ( step) ):
		### create a .txt with average consensus Qv for each sequence
		f = FastqReader(reseq_output+".fastq")
		with open(self.args.proj+"/reseq/quality_summary_%s.txt" % ( step), "w" ) as of:
		    for r in f:
			print>>of, "%s %s" % (r.name, sum([ (ord(x) - 33 ) for x in r.quality ])/float(len(r.quality)))

	    ### separate out resequenced HLA sequences from the nonspecific ones, by sequence name, store in two different files
	    ### TODO: read fastq output, parse consnesus QV, store somewhere the annotator can get at it
	    target_seq_names=[]
            for record in FastaReader(self.args.proj+"/phased/haplotype_consensus_sequences.fasta"):
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
		return 0
	    try:
                os.mkdir(self.args.proj+"/annotate")
            except:
                pass
	    self.logger.info("Annotation begun")
	    locus_dict={}	
	    nseqs=0
	    with open(self.args.proj+"/phased/phasr_output_seqs.txt", "r") as f:
		for line in f:
		    locus_dict[line.strip().split()[0]] = line.strip().split()[1]
		    nseqs+=1
	    self.logger.info("( %s ) sequences to annotate." % (nseqs) )
	    MSA_fn_dict={}; MSA_cDNA_fn_dict={}
	    MSA_info_fn_dict={}; MSA_cDNA_info_fn_dict={}
	    tmp_fn=self.args.proj+"/annotate/tmp.fasta"
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
	    f = FastqReader(self.args.proj+"/reseq/resequenced_hap_con.fastq")
	    if os.path.isfile(self.args.proj+"/annotate/gDNA_allele_calls.txt"):
		os.remove(self.args.proj+"/annotate/gDNA_allele_calls.txt")
	    if os.path.isfile(self.args.proj+"/annotate/resequenced_hap_con_cDNA.fasta"):
		os.remove(self.args.proj+"/annotate/resequenced_hap_con_cDNA.fasta")
	    if os.path.isfile(self.args.proj+"/annotate/cDNA_allele_calls.txt"):
		os.remove(self.args.proj+"/annotate/cDNA_allele_calls.txt")
	    for r in f:
		write_fasta([r], tmp_fn, "w") 
		name = r.name.split("|")[0]
		locus = locus_dict[name]
		output = self.args.proj+"/annotate/"+name+".afa"
		self.logger.info("Processing ( %s ) from locus ( %s )." % (name, locus) )

		### read in profile features
		MSA_info_fn = MSA_info_fn_dict[locus]
		MSA_info = {}
		info = namedtuple('info', 'canonical_pos, feature, codon' )
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
	    locus_counts = Counter(locus_dict.values())
	    loci_to_compare=[]
	    for item in locus_counts.iteritems():
		if item[1] == 2:
		    loci_to_compare.append(item[0])

	    phase_writer = GffWriter(self.args.proj+"/annotate/phase.gff")
	    phase_writer.writeMetaData('pacbio-variant-version', '1.4')

	    for locus in loci_to_compare:
		seqs_to_compare=[]
		for item in locus_dict.iteritems():
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
		
	    return 0
				
if __name__ == '__main__':
	if len(sys.argv) < 2:
		raise SystemExit("See `hla_pipeline.py --help` for usage information")
	HLA = HLApipeline()
        HLA.subread_extraction()
        HLA.phasing()	
	HLA.resequence()
	HLA.annotate()
