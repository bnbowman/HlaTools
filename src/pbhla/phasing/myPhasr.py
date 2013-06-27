#!/home/UNIXHOME/jquinn/HGAP_env/bin/python

#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$

import logging
import argparse
import sys
import os
import subprocess
import multiprocessing
import pkg_resources
import subprocess
import itertools
import signal
import shutil
import random
import string
import tempfile

from math import floor, log, ceil
from collections import namedtuple

from pbcore.io.FastaIO import FastaReader
from pbtools.pbdagcon.aligngraph import *
from pbtools.pbdagcon.utils import *

from myPhasrUtils import *

__p4revision__ = ""
__p4change__ = ""
revNum = int(1)
changeNum = int(0)
__version__ = "0.1.0" 

rmap = dict(zip("ACGTN-","TGCAN-"))
fastar = namedtuple('fastar', 'name, sequence')

class Feature(object):
    def __init__(self, name, tmp_dir):
	self.name = name
	self.metric = 0
	self.pctsimilarity = 0
	self.entropy = 0
	self.h1_reads = []
	self.h1_alns = []
	self.h1_con = ''
	self.h2_reads=[]
	self.h2_alns=[]
	self.h2_con=''
	self.h1_backbone=''
	self.h2_backbone=''
	self.mismatch = None
	self.insertion = None
	self.deletion = None
	self.aln_portion = None
	self.flag=0 ## feature is considered dead if flag = 1	
	self.tmp_dir = tmp_dir

    def __str__(self):
	output=[]	
	output.append("Name: %s" % (self.name) )
	output.append("Metric: %s" % (self.metric) )
	output.append("Entropy: %s" % (self.entropy) )
	output.append("PctID: %s" % (self.pctsimilarity) )
	if self.aln_portion:
	    output.append("AlnRatio: %s" % (self.aln_portion) )
	if self.insertion:
	    output.append("(MM: %s, ID: %s)" % (self.mismatch, (self.insertion+self.deletion) ) ) 
	output.append("Len: ( %s, %s )" % (len(self.h1_con), len(self.h2_con) ) )
	output.append("Size: ( %s, %s )" % (len(self.h1_alns), len(self.h2_alns) ) )
	if self.flag:
	    output.append("Flagged")
	string = " ".join(output)
	return string

    def refine(self, fasta_fn, n_iter):
	assert n_iter > 0
	for i in xrange(n_iter):	
	    if self.flag: continue
	    seqs_fn = self.write_seqs(self.tmp_dir, split=False)
	    tmp_fn = os.path.join( self.tmp_dir, (self.name+"_alignment.m5") )
	    subprocess.check_output("blasr -bestn 1 -m 5 %s %s -out %s" % (fasta_fn, seqs_fn, tmp_fn), shell=True )
	    self.h1_alns = get_aln_array( simple_align_hit_iterator(tmp_fn, "h1" ), max_num_reads=9999)	
	    self.h2_alns = get_aln_array( simple_align_hit_iterator(tmp_fn, "h2" ), max_num_reads=9999)
	    new_h1_consensus = make_template_from_alns( self.h1_alns, self.h1_con, combo_entropy = False)
	    new_h2_consensus = make_template_from_alns( self.h2_alns, self.h2_con, combo_entropy = False)
	    self.h1_backbone = self.h1_con; self.h2_backbone = self.h2_con
	    self.h1_con = new_h1_consensus.sequence
	    self.h2_con = new_h2_consensus.sequence
	    self.evaluate_pct_id()
	    self.entropy = (new_h1_consensus.numerator + new_h2_consensus.numerator )/float(new_h1_consensus.denominator + new_h2_consensus.denominator)
	    self.metric = (self.entropy*self.pctsimilarity)

	    ### cleanup
	    os.remove(tmp_fn)
	    os.remove(seqs_fn)
	return 0

    def finalize(self, fasta_fn, backbone, n_refinement):
	assert os.path.isfile(fasta_fn)
	if self.flag: return 0	
	seqs_fn = self.write_seqs(self.tmp_dir, split=False)
	tmp_fn = os.path.join( self.tmp_dir, (self.name+"_alignment.m5") )
	subprocess.check_output("blasr -bestn 1 -m 5 %s %s -out %s" % (fasta_fn, seqs_fn, tmp_fn), shell=True )
	self.h1_alns = get_aln_array( simple_align_hit_iterator(tmp_fn, "h1" ), max_num_reads=9999)
	self.h2_alns = get_aln_array( simple_align_hit_iterator(tmp_fn, "h2" ), max_num_reads=9999)	

	### write out the reads
	reads_h1_fn = os.path.join( self.tmp_dir, (self.name+"_h1_reads.fasta") )
	read_names = [ x[2].split("/")[0] for x in self.h1_alns]
	f=FastaReader(fasta_fn)
	with open( reads_h1_fn, "w") as of:
	    for r in f:
		if r.name in read_names: 
			print >>of, ">"+r.name
			print >>of, r.sequence
	reads_h2_fn = os.path.join( self.tmp_dir, (self.name+"_h2_reads.fasta") )
	read_names = [ x[2].split("/")[0] for x in self.h2_alns]
        f=FastaReader(fasta_fn)
        with open( reads_h2_fn, "w") as of:
            for r in f:
                if r.name in read_names:
                        print >>of, ">"+r.name
                        print >>of, r.sequence


	self.h1_backbone = backbone
	self.h1_con = get_good_consensus(self.h1_alns, backbone, reads_h1_fn, min_iteration = n_refinement)
	self.h2_backbone = backbone	
	self.h2_con = get_good_consensus(self.h2_alns, backbone, reads_h2_fn, min_iteration = n_refinement)
	self.evaluate_pct_id()

	### cleanup
	os.remove(tmp_fn)
	os.remove(seqs_fn)

	### return filenames of the reads
	return reads_h1_fn, reads_h2_fn

    def normalize(self):
	### by normalize i mean cut the sequences down only to the overlapping portion. This will prevent contig
	### effects for larger sequences, when reads from the same side cluster instead of from the same allele
	if self.flag: return 0
        if self.h1_con == self.h2_con:
            self.flag = 1
            return 0
	seq_fns = self.write_seqs(self.tmp_dir, split=True)
	try:
	    alignment = parse_blasr( subprocess.check_output("blasr -bestn 1 %s %s -m 5" % (seq_fns[0], seq_fns[1]), shell=True), mode=5)
	    alignment = alignment[0]
	    self.h1_con = alignment.qseq.replace('-', '')	
	    self.h2_con = alignment.tseq.replace('-', '')
	except:
            self.flag=1
	    for fn in seq_fns:
		os.remove(fn)
	    return 0
        for fn in seq_fns:
            os.remove(fn)	

    def evaluate_pct_id(self):
	if self.flag: return 0
	if self.h1_con == self.h2_con: 
	    self.flag = 1
	    return 0
	seq_fns = self.write_seqs(self.tmp_dir, split=True)
	try:
	    alignment = parse_blasr( subprocess.check_output("blasr -bestn 1 %s %s -m 4" % (seq_fns[0], seq_fns[1]), shell=True), mode=4)
	    self.pctsimilarity = float(alignment[0].pctsimilarity)
	    alignment = parse_blasr( subprocess.check_output("blasr -bestn 1 %s %s -m 5" % (seq_fns[0], seq_fns[1]), shell=True), mode=5)
	    alignment = alignment[0]
	    ### the percentage of mismatches due to ins, del, mismatch etc
	    self.mismatch = float(alignment.nmis)/len(alignment.matchvector)
	    self.insertion = float(alignment.nins)/len(alignment.matchvector)
	    self.deletion = float(alignment.ndel)/len(alignment.matchvector)
	    ### how long is the alignment compared to how long it *could* be
	    self.aln_portion = ( len(alignment.matchvector) / float(min([ len(self.h1_con), len(self.h2_con) ])) )
	except:
	    self.flag=1
	for fn in seq_fns:
	    os.remove(fn)
	
    def write_seqs(self, outdir, split = True):
	assert os.path.isdir(outdir)
	if split:
	    with open(os.path.join(outdir, (self.name+"_h1_con.fasta")), "w") as of:
		print >>of, ">h1"
		print >>of, self.h1_con
	    with open(os.path.join(outdir, (self.name+"_h2_con.fasta")), "w") as of:
		print >>of, ">h2"
		print >>of, self.h2_con
	    return ( os.path.join(outdir, (self.name+"_h1_con.fasta") ), os.path.join(outdir, (self.name+"_h2_con.fasta") ) )
	elif not split:
	    with open(os.path.join(outdir, (self.name+"_h_con.fasta")), "w") as of:
		print >>of, ">h1"
                print >>of, self.h1_con
		print >>of, ">h2"
                print >>of, self.h2_con
	    return os.path.join(outdir, (self.name+"_h_con.fasta"))

class Phasr(object):
    """
    Tool for separating out different alleles 
    """

    def __init__(self, inputFile=None, refFile=None):
        if inputFile is None or refFile is None:
            self.initializeFromArgs()
        else:
            pass
        self.initializeLogger()
        self.finishInitialization()

    def initializeFromArgs(self):
        parser = argparse.ArgumentParser(description = "Phase Pac Bio CLRs.")

        add = parser.add_argument
	add('fasta_fn', metavar='input.fasta', help='an input fasta file')
	add('--ref_fn',  metavar='ref.fasta', help='a reference fasta file')
	add('--output', metavar='output.fasta', dest='output_fn', 
	                    default=os.path.join( os.getcwd(), "h_consensus.fasta"), 
	                    help='Consensus output filename.')
	add('--sample_size', default=6, metavar='6', 
	                    type=int, dest='sample_size', 
	                    help='The maximum number of reads used for consensus')
	add('--sample_number', type=int, default=6,
	                    dest='sample_number', metavar='6',
	                    help='The maximum number of reads used for consensus')
	add('--min_cluster_size', type=int, default=10, 
	                    dest='min_cluster_size', metavar='10', 
	                    help='The minimum number of reads needed to define a cluster.')
	add('--min_cluster_divergence', type=float, default=99.5, 
	                    dest='min_cluster_divergence', metavar='99.5', 
	                    help='The minimum divergence between two ' + \
	                    'clusters for them to be separated.')
	add('--max_recursion_level', type=int, default=2, 
	                    dest='max_recursion_level', metavar='2', 
	                    help='Once recursion reaches this level, ' + \
	                    'consensus will be made from the current grouping ' + \
	                    'and terminate. If set to 0, will just make' + \
	                    'consesus from multifasta.')
	add('--n_refinement', type=int, default=2, 
	                    dest='n_refinement', metavar='2', 
	                    help='Number of dagcon iterations for building ' + \
	                    'consensus. Making this number higher will cause ' + \
	                    'slower execution. Max useful setting is around ' + \
	                    '~5, will help deal with structural rearrangements.')
	add('--fullpass', action='store_true', dest='fullpass', 
	                    help='Only fullpass reads allowed. Will look for' + \
	                    ' the "fp" tag in fasta sequence names.')
	add('--score_floor', type=float, default=0.0, 
	                    dest = 'score_floor', metavar='0.0', 
	                    help='Min percent id for construction of Max Divergent Features.')
	add('--log', action='store_true', dest='log', 
	                    help = 'Create log file.')
	add('--max_input_reads', type=int, default=400, 
	                    dest='max_input_reads', metavar='400', 
	                    help='Maximum number of reads to use')
	add('--min_read_length', type=int, default=500, 
	                    dest='min_read_length', metavar='500',
	                    help = 'Minimum read size' )
	add('--max_num_proc', type=int, default=4, 
	                    dest='max_num_proc', metavar='4',
	                    help='Maximum number of subprocesses to spawn. ' + \
	                    'Total processes will be this number + 1 for ' + \
	                    'the parent process.')
	self.args = parser.parse_args()
        self.output_dir = os.path.dirname(self.args.output_fn)

    def initializeLogger(self):
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)
	if self.args.log:
	    logFile = os.path.join( self.output_dir, "phasr.log" )
	    h1 = logging.FileHandler( logFile )
	elif not self.args.log:
	    h1 = logging.FileHandler("/dev/null")
	f = logging.Formatter("%(processName)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
	h1.setFormatter(f) 
	h1.setLevel(logging.INFO) 
	logger.addHandler(h1)
	h2 = logging.StreamHandler(stream=sys.stdout)
	h2.setFormatter(f)
	h2.setLevel(logging.INFO)
	logger.addHandler(h2)
	self.logger = logger
	self.logger.info("phasr envoked: %s" % " ".join(sys.argv) )
	if self.args.log:
	    self.logger.info("Logging in %s" % (os.path.join( os.path.dirname(self.args.output_fn), "hla_pipeline.log" )) )


    def finishInitialization(self):
	self.fasta_stack = [] ### list of fasta files to process
	self.consensus_dictionary={}
	self.hap_cons=[]
        ### catch signals
	signal.signal(signal.SIGINT, self.signal_handler)	
	### check arguments
	assert os.path.isfile(self.args.fasta_fn)	
	if self.args.ref_fn is not None:
	    assert os.path.isfile(self.args.ref_fn)
	    self.ref_fn = os.path.abspath(self.args.ref_fn)
	assert os.path.exists(os.path.dirname( self.args.output_fn ) )
	if os.path.isfile(self.args.output_fn):
	    self.logger.info("%s already exists, overwriting.." %  (self.args.output_fn) )	
        self.tmp_dir = tempfile.mkdtemp()
        self.logger.info("Tmp Dir is %s" % self.tmp_dir )
        self.prepareInputFasta()

    def prepareInputFasta(self):
	### prepare input fasta
	self.input_fn = os.path.join( self.tmp_dir, os.path.basename(self.args.fasta_fn) )
	total_reads=0 
	used_reads=0
        f = FastaReader(self.args.fasta_fn) 
        with open(self.input_fn, "w") as output:
            for record in f:
                total_reads+=1
                if self.args.fullpass:
                    if 'fp' not in record.name: continue
                if self.args.min_read_length:
                    if len(record.sequence) < self.args.min_read_length: continue
                used_reads+=1
                print >>output, '>' + record.name
                print >>output, record.sequence
	self.logger.info("%s/%s reads remain after filtering." % (used_reads, total_reads) )
	### limit the number of reads used in the algorithm if specified
	if used_reads > self.args.max_input_reads:
	    f = FastaReader(self.input_fn)
	    fp_input_reads=[]; non_fp_input_reads=[];
	    for r in f:
		if r.name.startswith('fp'):
		    fp_input_reads.append(r)
		else:
		    non_fp_input_reads.append(r)
	    if (len(non_fp_input_reads)+len(fp_input_reads)) > int(self.args.max_input_reads):
		if len(fp_input_reads) < int(self.args.max_input_reads):
		    non_fp_reads_size = (int(self.args.max_input_reads)-len(fp_input_reads))
		    non_fp_input_reads = random.sample(non_fp_input_reads, non_fp_reads_size)
		    fp_input_reads.extend(non_fp_input_reads)
		    input_reads = fp_input_reads
		else:
		    input_reads = random.sample(fp_input_reads, int(self.args.max_input_reads))
		with open(tmp_fn, "w") as of:
		    for r in input_reads:
			print >>of, ">"+r.name
			print >>of, r.sequence
	    os.remove(self.input_fn); shutil.move(tmp_fn, self.input_fn)
	    self.logger.info("%s reads will be used." % (self.args.max_input_reads)  )
	### initialize the stack with the input fasta file
	self.fasta_stack.append((self.input_fn, 0))
	print self.fasta_stack

    def signal_handler(signal, frame):
        self.logger.info("Recieved SIGINT.")
        self.cleanup()
        raise SystemExit
 
    ##TODO
    def cleanup(self):
	shutil.rmtree(self.tmp_dir)

    def getVersion(self):
        return __version__

    def generate_haplotype_consensus(self): 
	    
	### define filenames
	tmp_fasta = os.path.join( self.tmp_dir, make_rand_string()+".fasta" )	
	tmp_file = os.path.join( self.tmp_dir, make_rand_string())	

	### get fasta and ref
	try:
	    input_fn, rec_level = self.fasta_stack.pop()
	except:
	    return 0

	self.logger.info("Now processing ( %s ) at level ( %s )." % (input_fn, rec_level) )
	f = FastaReader(self.ref_fn)
	for r in f: backboneSeq = r.sequence; break

	### normalize fasta and get alns

	normalize_fasta(input_fn, self.ref_fn, tmp_fasta)
	os.system("blasr %s %s -m 5 -out %s" % (tmp_fasta, self.ref_fn, tmp_file))
	alns = get_aln_array( simple_align_hit_iterator(tmp_file), max_num_reads=9999)
	shutil.move(tmp_fasta, input_fn)
	os.remove(tmp_file)

	try:
	    ### see if we already made a consensus from this subset
	    consensus = self.consensus_dictionary[os.path.abspath(input_fn)]
	except KeyError:
	    ### or make initial consensus using all reads from this subset
	    self.logger.info("%s: Creating initial consensus" % (rec_level) )
	    consensus = get_good_consensus(alns, backboneSeq, input_fn)
	    self.consensus_dictionary[os.path.abspath(input_fn)] = consensus 
	init_seq_length = len(consensus)
	self.logger.info("%s: Initial sequence is of length ( %s )" % (rec_level, init_seq_length) )

	### initial conditions under which phasing will cease and the initial consensus will be returned
	if len(alns) < self.args.min_cluster_size or rec_level >= self.args.max_recursion_level:
	    self.logger.info("%s: Unable to further split this group. Will return initial consensus." % ( rec_level ) )
	    self.hap_cons.append(fastar._make( [ make_rand_string(), consensus ] ))
	    self.logger.info("%s: %s" % (rec_level, consensus))
	    return 0

	self.logger.info("%s: Phasing ( %s ) alns." % ( rec_level, len(alns) ) )

	### call workers to come up with Max Divergent Features using different
	### random samples of these alignments
	### throw the multiprocess objects in a dict then run them
	process_dict={}
	manager = multiprocessing.Manager()
	feature_list = manager.dict()

	for k in xrange(self.args.sample_number):	
	    while process_status(process_dict) >= self.args.max_num_proc:
		pass
	    process_dict[k] = multiprocessing.Process(target=create_feature, args=( alns, backboneSeq, feature_list, self.args.sample_size, init_seq_length, self.args.score_floor, self.tmp_dir, input_fn, self.args.n_refinement ))
	    process_dict[k].start()
	for item in process_dict.itervalues():
            item.join()

	feature_list = dict(feature_list)
	def ranking_function(feature):
	    aln_sizes = sorted([len(feature.h1_alns), len(feature.h2_alns)])
	    return float(aln_sizes[1])/float(aln_sizes[0])
	feature_ranking = sorted([ item for item in feature_list.itervalues() ], key = ranking_function )
	successful_features = []

	unflagged_features=0
	for feature in feature_ranking:
	    self.logger.info("%s: Created Feature: %s" % (rec_level, feature) )
	    if feature.flag: continue
	    unflagged_features+=1

	### exit condition 
	if not unflagged_features:
	    self.logger.info("%s: Unable to further split this group. Will return initial consensus." % ( rec_level ) )
	    self.hap_cons.append(fastar._make( [ make_rand_string(), consensus ] ))
	    self.logger.info("%s: %s" % (rec_level, consensus))
	    return 0 

	### after iterating though the features and comparing their metrics, the one with the lowest metric is chosen
	### we refine it further, and then impose a few conditions on it
	### if it meets all of these stricter conditions, the function will exit and return as output two fasta files
	### which are subsets of the input fasta file
	while True:
	    try:
		current_feature=feature_ranking.pop(0)
	    except IndexError:
		self.logger.info("%s: Unable to further split this group. Will return initial consensus." % ( rec_level ) )
		self.hap_cons.append(fastar._make( [ make_rand_string(), consensus ] ))
		self.logger.info("%s: %s" % (rec_level, consensus))
		return 0
	    if current_feature.flag: continue
	    if len(current_feature.h1_alns) <= self.args.min_cluster_size or len(current_feature.h2_alns) <= self.args.min_cluster_size:
		continue
	    self.logger.info("%s: Processing Feature: %s" % (rec_level, current_feature) )
	    ### align all reads back to the feature and generate iterative dagcon consensus starting from backbone
	    read_subset1_fn, read_subset2_fn = current_feature.finalize(input_fn, backboneSeq, self.args.n_refinement)
	    ### sanity check the clustering results
	    if len(current_feature.h1_alns) <= self.args.min_cluster_size or len(current_feature.h2_alns) <= self.args.min_cluster_size:
		self.logger.info("%s: Feature ( %s ) failed due to small cluster size <= ( %s )." % ( rec_level, current_feature.name, self.args.min_cluster_size) )
		continue
	    if current_feature.pctsimilarity >= self.args.min_cluster_divergence:
		self.logger.info("%s: Feature ( %s ) failed due to high consensus sequence identity." % ( rec_level, current_feature.name) )
		continue
	    self.logger.info("%s: Feature ( %s ) will be used, recursion will descend." % ( rec_level, current_feature.name) )
	    self.logger.info("%s: %s" % (rec_level, current_feature) )
	    self.logger.info("%s: Sequence 1 of the successful feature is: %s" % ( rec_level, current_feature.h1_con ) )
	    self.logger.info("%s: Sequence 2 of the successful feature is: %s" % ( rec_level, current_feature.h2_con ) )

	    ### put the subgroups on the stack
	    self.fasta_stack.insert(0, (read_subset1_fn, rec_level+1))
	    self.fasta_stack.insert(0, (read_subset2_fn, rec_level+1))
	    ### save generating consensus from the same group of reads more than once
	    self.consensus_dictionary[os.path.abspath(read_subset1_fn)] = current_feature.h1_con
	    self.consensus_dictionary[os.path.abspath(read_subset2_fn)] = current_feature.h2_con
	    return 0

    def run(self):
	if self.args.ref_fn == None: ### if denovo mode is chosen a read is used as the template
	    try:
		ref=fastar._make(best_template_by_blasr(input_fasta_name))
	    except:
		self.cleanup()
		self.logger.info("De novo template selection failed. Exiting..")
		return 0
	    write_fasta(ref, os.path.join(self.tmp_dir, "btbb.fasta") )
	    self.args.ref_fn = os.path.join(self.tmp_dir, "btbb.fasta")

	while 1:
	    self.generate_haplotype_consensus()
	    if len(self.fasta_stack) == 0: break
	if len(self.hap_cons) > 0: write_fasta( self.hap_cons, self.args.output_fn)
	self.logger.info("( %s ) sequences output to ( %s )" % ( len(self.hap_cons), self.args.output_fn ) )
	seq_to_read_fn = dict([[v,k] for k,v in self.consensus_dictionary.items()]) ##TODO: hash sequence 
	for seq in self.hap_cons:
	    read_fn = seq_to_read_fn[seq.sequence]
	    shutil.copyfile(read_fn, os.path.join(os.path.dirname(self.args.output_fn), seq.name+".fasta"))	

	self.logger.info("Process complete")
	self.cleanup()

def create_feature(alns, backboneSeq, feature_list, sample_size, init_seq_length, score_floor, tmp_dir, input_fn, n_refinement):
    while 1: ### give a collision-impossible name to this feature
	rands=make_rand_string() 
	temp1_fn =os.path.join(tmp_dir, rands+"_h1.fasta")
	temp2_fn = os.path.join(tmp_dir, rands+"_h2.fasta")
	if rands not in feature_list: break
    reads=random.sample(alns, sample_size)
    worstscore = float(10000000000.0) ### aim to minimize this number
    output = None
    bestentropy = float(999) ### aim to minimize this number
    worstidentity=float(100)
    ### begin iterating through all possible groupings of this subset of reads
    for i in range(1,int((floor(len(reads)/2))+1)):
	for combo in itertools.combinations(reads, i):
	    icombo = list(set(reads).difference(combo))
	    ### get the two consensus templates
	    templates = []
	    for subset in [combo, icombo]:
		template = make_template_from_alns(subset, backboneSeq, max_num_reads = 10)
		templates.append(template)
	    try:
		entropy = sum([ (lambda x: x.numerator)(x) for x in templates ])/float(sum([ (lambda x: x.denominator)(x) for x in templates ]))
	    except:
		entropy = 1
	    ### we need to write the consensus sequences to file in order to run blasr, and get percent ID
	    ### TODO: find away around doing this, its slow!
	    with open(temp1_fn, "w") as of1:
		print >>of1, ">h1"
		print >>of1, templates[0].sequence
	    with open(temp2_fn, "w") as of2:
		print >>of2, ">h2"
		print >>of2, templates[1].sequence 
	    ### get percent ID 
	    try:
		blasr_output = parse_blasr(subprocess.check_output("blasr -bestn 1 -m 4 %s %s" % ( temp1_fn, temp2_fn), shell=True), 4)
		score = float(blasr_output[0].pctsimilarity)
	    except:
		continue
	    ### set a lower bound on similarity, we dont want to return garbage
	    if score <= score_floor:
		continue
	    ### calculate our metric (percent ID * entropy), if the metric is lower than any previous combination
	    ### write it to memory
	    if (score * entropy) < worstscore:
		output = combo
		worstscore = (score * entropy)
		bestentropy = entropy
		worstidentity = score
		best_template_pair = templates
    ### return the feature to the process manager
    if output == None:
	return 0
    ioutput = list(set(reads).difference(output))
    created_feature = Feature(rands, tmp_dir)
    created_feature.metric = worstscore; created_feature.entropy = bestentropy
    created_feature.h1_con = best_template_pair[0].sequence; created_feature.h2_con = best_template_pair[1].sequence 
    created_feature.h1_alns = output
    created_feature.h2_alns = ioutput
    created_feature.h2_backbone = backboneSeq
    created_feature.h1_backbone = backboneSeq
    created_feature.evaluate_pct_id()

    ### now align all reads to feature and refine 
    created_feature.refine( input_fn, n_refinement )
    created_feature.normalize()
	
    feature_list[rands] = created_feature 

    ### cleanup
    os.remove(temp1_fn); os.remove(temp2_fn)

if __name__ == '__main__':    
    sys.exit(Phasr().run())
