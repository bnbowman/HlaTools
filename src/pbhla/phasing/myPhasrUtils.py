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

from math import floor, log, ceil
from collections import namedtuple

from pbcore.io.FastaIO import FastaReader
from pbtools.pbdagcon.aligngraph import *
from pbtools.pbdagcon.utils import *

__p4revision__ = ""
__p4change__ = ""
revNum = int(1)
changeNum = int(0)
__version__ = "0.1.0" 

rmap = dict(zip("ACGTN-","TGCAN-"))
fastar = namedtuple('fastar', 'name, sequence')

def write_fasta(fasta_obj, outfile, mode = "w"):
    if isinstance(fasta_obj, list):
	with open(outfile, mode) as of:
	    for r in fasta_obj:
		if r.name[0] != ">":
		    print >>of, ">"+r.name
		elif r.name[0] == ">":
		    print >>of, r.name
		print >>of, r.sequence
    else:
	r = fasta_obj
	with open(outfile, mode) as of:
	    if r.name[0] != ">":
		print >>of, ">"+r.name
	    elif r.name[0] == ">":
		print >>of, r.name
	    print >>of, r.sequence

def fasta_size(fasta):
    try:
        f = FastaReader(fasta)
        count=0
        for r in f:
            count+=1
        return count
    except:
        return None

def extract_sequence(fasta, names):
    f = FastaReader(fasta)
    if isinstance(names, str):
        for r in f:
            if r.name == names:
                return r
    elif isinstance(names, list):
        output=[]
        for r in f:
            if r.name in names:
                output.append(r)
        return output

def process_status(process_dict):
        alive_processes=0
        for p in process_dict.itervalues():
            if p.is_alive():
                alive_processes+=1
        return alive_processes

def parse_blasr(output, mode, strip_query_names = True):
    parsed_output=[]
    if mode == 1:
        entry = namedtuple('entry', 'qname, tname, qstrand, tstrand, score, pctsimilarity, tstart, tend, tlength, qstart, qend, qlength, ncells')
    if mode == 4:
        entry = namedtuple('entry', 'qname, tname, score, pctsimilarity, qstrand, qstart, qend, qseqlength, tstrand, tstart, tend, tseqlength, mapqv, ncells, clusterScore, probscore, numSigClusters')
    if mode == 5:
        entry = namedtuple('entry', 'qname, qlength, z1, qalength, qstrand, tname, tlength, z2, talength, tstrand, score, nmatch, nmis, nins, ndel, zscore, qseq, matchvector, tseq')

    output = output.strip().split("\n")
    output = [ x.split() for x in output ]
    if strip_query_names:
	new_output=[]
	for line in output:
	    line[0] = line[0].split("/")[0]
	    new_output.append(line)
	output = new_output

    for line in output:
	alignment = entry._make(line)
        parsed_output.append(alignment)
    return parsed_output

def make_rand_string(minlength=6,maxlength=8):
    length=random.randint(minlength,maxlength)
    letters=string.ascii_letters+string.digits # alphanumeric, upper and lowercase
    return ''.join([random.choice(letters) for _ in range(length)])

def best_template_by_blasr(fasta_fn, len_threshold = 300, min_number_reads = 1, rank_reads = False):
    f = FastaReader(fasta_fn)
    read_dict = {}
    for r in f:
        r_id = r.name.split("/")[0]
        read_dict[r_id] = r.sequence

    rtn = os.system("blasr -bestn 20 -nCandidates 100 -m 1 %s %s -out %s" % (fasta_fn, fasta_fn, fasta_fn+".saln"))
    if rtn != 0:
        return None

    scores = {}
    with open(fasta_fn + ".saln") as f:
        for l in f:
            l = l.strip().split()
            l = l[:9]
            r1, r2 = l[:2]
            r1 = r1.split("/")[0]
            r2 = r2.split("/")[0]
            if r1 == r2:
                continue
            if int(l[7]) - int(l[6]) < len_threshold : continue
            scores.setdefault(r1,[])
            scores[r1].append(int(l[4]))
            #scores[r1].append(-float(l[5]))
    score_array = []
    for r, s in scores.items():
        score_array.append( (np.mean(s), r) )
    score_array.sort()
    if rank_reads:
	return score_array

    if min_number_reads != None:
        if len(score_array) < min_number_reads:
            raise AlignGraphUtilError("not enough number of reads in best_template_by_blasr")
        else:
            return score_array[0][1], read_dict[score_array[0][1]]
    else:
        return r_id, read_dict[r_id]

def construct_aln_graph_from_aln_array(alns, 
			backboneSeq,
			max_num_reads = None,
			remove_in_del = False,
			ref_group = None,
			min_length = None):

    g = AlnGraph(backboneSeq)

    i = 0
    for aln in alns:
	rId = aln[2]
	aln = aln[0:2]
	if i <= max_num_reads:
	    try:
		g.add_alignment( aln, "%s" % rId)
		i+=1
	    except:
		pass
	else:
	    break	
    return g

def normalize_fasta(fasta_file, ref_file, out_file):
    f = FastaReader(fasta_file)
    with open(out_file, "w") as of:
        for r in f:
            r_id = "%s" %  r.name
            print >>of, ">"+r_id
            seq = r.sequence.upper()
            print >>of, seq 

    output = subprocess.check_output("blasr -bestn 1 -m 1 %s %s" % ( out_file, ref_file ), shell=True)
    direction = {}
    output = output.strip().split("\n")
    for l in output:
        l = l.strip().split()
        rId = l[0].split("/")[0]
        if l[2] != l[3]:
            direction[rId] = "-"
        else:
            direction[rId] = "+"

    f = FastaReader(out_file)
    outData = []
    for r in f:
        r_id = "%s" % r.name
        outData.append(">"+r_id)
        seq = r.sequence.upper()
        if direction != None:
            if direction.get(r_id, "+") != "+":
                seq = "".join([rmap[c] for c in seq[::-1]])
        outData.append(seq)
    with open(out_file,"w") as of:
        print >>of, "\n".join(outData)

def get_good_consensus(alns, backboneSeq, read_fn,
                  hp_correction = True,
                  min_iteration = 2,
                  max_num_reads = 1000,
                  entropy_th = 0.65,
		  min_cov = 0,
		  consensus_seq_name = 'consensus'):
    while 1:
	tmp_fn = os.path.join( os.path.dirname(read_fn), (make_rand_string()+"_con.fasta") )
	if not os.path.isfile(tmp_fn) : break
    g = construct_aln_graph_from_aln_array(alns, backboneSeq, max_num_reads = max_num_reads, remove_in_del = False)
    s,c = g.generate_consensus(min_cov = min_cov)
    output = [ fastar._make([ consensus_seq_name,  s.upper()] ) ]

    assert os.path.isdir(os.path.dirname(tmp_fn))
    write_fasta( output, tmp_fn )	

    for j in range(2):
	for i in range(min_iteration-2):
	    g = constructe_aln_graph_from_fasta(read_fn, tmp_fn, max_num_reads = max_num_reads, remove_in_del = False)
	    s,c = g.generate_consensus(min_cov = min_cov)
	    output = [ fastar._make([consensus_seq_name,  s.upper()]) ]
	    write_fasta( output, tmp_fn)	

	if hp_correction:
	    g = constructe_aln_graph_from_fasta(read_fn, tmp_fn, max_num_reads = max_num_reads, remove_in_del = False)
	    s = detect_missing(g, entropy_th = entropy_th)
	    output = [ fastar._make([consensus_seq_name,  s.upper()]) ]
	    write_fasta( output, tmp_fn)

	g = constructe_aln_graph_from_fasta(read_fn, tmp_fn, max_num_reads = max_num_reads, remove_in_del = False)
	s,c = g.generate_consensus(min_cov = min_cov)
	s = mark_lower_case_base(g, entropy_th = entropy_th)
	output = [ fastar._make([consensus_seq_name,  s.upper()]) ]
	write_fasta( output, tmp_fn)
    os.remove(tmp_fn)
    return s

def make_template_from_alns(aln_array, backboneSeq, 
                  hp_correction = True,
                  min_iteration = 1, 
                  max_num_reads = 150,
                  entropy_th = 0.65,
		combo_entropy=True):
    template_info = namedtuple("template_info", "sequence, numerator, denominator") # this function can produce two types of 'entropy' measures
    # both of which have a numerator and a denominator
    g = construct_aln_graph_from_aln_array(aln_array, backboneSeq, max_num_reads = max_num_reads)
    seq_length=len(backboneSeq)
    s,c = g.generate_consensus(min_cov = 0)
    nodes = sorted_nodes(g)
    denominator=0; numerator=0
    if combo_entropy:
	### if we have a small number of reads we will calculate entropy
	### as the proportion of nodes in the graph which have 100% of the reads in the set
	### passing through them
	for n in nodes:
	    if len(aln_array) > 1: ### special case of a set of one read, consider the backbone to be a read
		if len(n.info) == 0:
		    continue
	    edge_list=[]
	    denominator+=1
	    for item in n._in_edges: edge_list.append(int(item.count))
	    if sum(edge_list) < len(aln_array): numerator+=1
    elif not combo_entropy:
	### if we are claculating entropy for a larger group of reads we do it differently
	### since it is likely that 0 nodes will have 100% of all reads going through them
	### the main idea here is to get a statistical estimate of the cross entropy
	### that will be relatively independent of cluster size
	### thus for each node we caclulate 1/cluster_size * -log( cov/max_cov)
	for n in nodes:
	    edge_list=[]
	    if len(aln_array) > 1: ### special case of a set of one read, consider the backbone to be a read
		if len(n.info) == 0:
		    continue
	    for item in n._in_edges: edge_list.append(int(item.count))
	    if sum(edge_list) == 0: continue
	    p = log(float(sum(edge_list)/float(len(aln_array)+1)))
	    q = (float(1/float(len(aln_array)+1)))
	    numerator+=(-1*(q*p))
	denominator=(-1 * log(1/float(len(aln_array)+1)) * (1/float(len(aln_array)+1)) * float(seq_length) * float(len(aln_array)+1))
    output = template_info._make( [s, numerator, denominator])
    return output
