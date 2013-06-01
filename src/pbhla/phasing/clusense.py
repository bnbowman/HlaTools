#! /usr/bin/env python

import os, sys, glob, logging
import pkg_resources
from math import log

import numpy as np
from pbcore.io import FastaReader
from pbtools.pbdagcon.c_aligngraph import *
from pbtools.pbdagcon.c_utils import construct_aln_graph_from_fasta
from pbtools.pbdagcon.vis_utils import output_read_network
from pbtools.pbdagcon.c_utils import sorted_nodes
from pbtools.pbdagcon.vis_utils import dump_sorted_node_data
from pbtools.pbdagcon.c_utils import best_template_by_blasr
from pbtools.pbdagcon.c_utils import clustering_read
from pbtools.pbdagcon.c_utils import get_subset_reads
from pbtools.pbdagcon.c_utils import read_node_vector
from pbtools.pbdagcon.c_utils import detect_missing

# Default values
MIN_GROUP = 25
THRESHOLD = 0.1
PREFIX = 'Unknown'

## prepare the reads
#!blasr test_HLA_A.fa u0038.fa -bestn 1 -out tmp.sam -noSplitSubreads -sam -nproc 32
#!cat tmp.sam | awk 'NF > 15 && length($10) > 2500 && i < 1000 {print ">"$1"\n"$10;i+=2}' > test.fa

def calculate_entropy( threshold ):
    return -1 * threshold * log(threshold) - (1-threshold) * log(1-threshold)

def get_consensus(read_fn, init_ref, consensus_fn, consens_seq_name, 
                  ENTROPY_TH,
                  hp_correction = False, 
                  min_iteration = 4, 
                  max_num_reads = 150,
                  entropy_th = 0.65,
                  min_cov = 8,
                  max_cov = 60,
                  nproc = 4,
                  mark_lower_case = False,
                  use_read_id = False):

    g = construct_aln_graph_from_fasta(read_fn, init_ref, max_num_reads = max_num_reads, remove_in_del = False, max_cov = max_cov, nproc = nproc)
    s, c = g.generate_consensus(min_cov = min_cov)


    with open(consensus_fn,"w") as f:
        print >>f, ">"+consens_seq_name
        print >>f, s.upper()

    if min_iteration > 1:
        for i in range(min_iteration-2):
            if len(s) < 100:
                return s
            g = construct_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = max_num_reads, remove_in_del = False, max_cov = max_cov, nproc = nproc, use_read_id = use_read_id)
            s, c = g.generate_consensus(min_cov = min_cov)
            with open(consensus_fn,"w") as f:
                print >>f, ">"+consens_seq_name
                print >>f, s.upper()

        if hp_correction:
            if len(s) < 100:
                return s
            g = construct_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = max_num_reads, remove_in_del = False, max_cov = max_cov, nproc = nproc)
            s = detect_missing(g, entropy_th = ENTROPY_TH)
            with open(consensus_fn,"w") as f:
                print >>f, ">"+consens_seq_name
                print >>f, s.upper()

        if len(s) < 100:
            return s
        g = construct_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = max_num_reads, remove_in_del = True, max_cov = max_cov, nproc = nproc)
        s, c = g.generate_consensus(min_cov = min_cov)
        if mark_lower_case:
            s = mark_lower_case_base(g, entropy_th = ENTROPY_TH)
        with open(consensus_fn,"w") as f:
            print >>f, ">"+consens_seq_name
            print >>f, s
    return g, s

def fetch_read(in_file, out_file, id_set):
    in_f = FastaReader(in_file)
    with open(out_file,"w") as out_f:
        for r in in_f:
            if r.name not in id_set:
                continue
            print >> out_f, ">"+r.name
            print >> out_f, r.sequence
    in_f.file.close()

def get_hp_run_counts(n):
    rtn = [0, 0]
    count = 0
    p_node = n.get_best_in_node()
    while p_node != None and  p_node.get_base() == n.get_base():
        count += 1
        n = p_node
        p_node = n.get_best_in_node()
        if count > 5:
            break
    rtn[0] = count
    count = 0
    n_node = n.get_best_out_node()
    while n_node != None and n_node.get_base() == n.get_base():
        count += 1
        n = n_node
        n_node = n.get_best_out_node()
        if count > 5:
            break
    rtn[1] = count
    return rtn

def is_hp_node(n):
    c1, c2 = get_hp_run_counts(n)
    if c1 + c2 >= 3:
        return True
    else:
        return False

def read_node_vector(g, ENTROPY_TH, entropy_th = 0.65):
    
    ne, hne = g.get_high_entropy_nodes(coverage_th = 0, entropy_th = ENTROPY_TH)
    node_to_entropy = dict( [ (v[1],v[2]) for v in ne ] ) 
    read_ids = set()
    
    for n in g.get_nodes().values():
        for r in n.get_info():
            read_ids.add(r[0]) #r[0]: read id, r[1]: read base
    
    backbone_node_to_reads = {}
    for n in g.get_nodes().values():
        bn = n.get_backbone_node()
        backbone_node_to_reads.setdefault(bn, set() )
        for r_id, r_pos in n.get_info():
            backbone_node_to_reads[bn].add( r_id )
    
    backbone_node_to_pos = g.get_backbone_node_to_pos()
    sn = g.get_sorted_nodes()
    high_entropy_nodes = [ (n, backbone_node_to_pos[n.get_backbone_node()], node_to_entropy[n]) \
                            for n in sn if n in node_to_entropy and node_to_entropy[n] > entropy_th]
    high_entropy_nodes = [ n for n in high_entropy_nodes if not is_hp_node(n[0]) ] 

    read_to_nodes = {}
    for r_id in read_ids:
        read_to_nodes[r_id] = [" "] * len(high_entropy_nodes)
        for i,n in enumerate(high_entropy_nodes):
            if r_id in [x[0] for x in n[0].get_info()]:
                read_to_nodes[r_id][i] = n[0].get_base()
            else:
                bs,be = g.get_read_range()[r_id]
                bp = backbone_node_to_pos[n[0].get_backbone_node()]
                if bp >= bs and bp < be:
                    read_to_nodes[r_id][i] = "*"
    return read_to_nodes, high_entropy_nodes

def partition_score(m, pos, index):
    m0 = m[m[:,pos] == 1 ,:]
    m1 = m[m[:,pos] == -1 ,:]
    m0 = m0[:,index]
    m1 = m1[:,index]
    return np.sum(np.sum(m0,0) *  np.sum(m1,0))

def col_entropy(m):
    sp = np.sum(m == 1,0) + 1
    sn = np.sum(m == -1, 0) + 1
    #sz = np.sum(m == 0, 0) + 1
    all_n = sp + sn
    # all_n = sp + sn + sz
    p1 = 1.0 * sp  / all_n
    p2 = 1.0 * sn  / all_n
    print len(p1), len(p2)
    #p3 = 1.0 * sz  / all_n
    #return -p1 * np.log(p1) - p2 * np.log(p2) - p3 * log(p3)
    return -p1 * np.log(p1) - p2 * np.log(p2)

def regroup(g1, g2, index):
    smap = {"A":1,"C":1,"G":1,"T":1,"*":-1," ":0}
    read_vector = g1 + g2
    m = {}
    m1 = np.zeros( (len(g1), len(index) ) )
    m2 = np.zeros( (len(g2), len(index) ) )
    g1_id = [x[0] for x in g1] 
    g2_id = [x[0] for x in g2] 
    all_IDs = g1_id + g2_id
    i1 = 0
    i2 = 0

    for ID, node in read_vector:
        vec = 1.0 * np.array([smap[c] for c in node])
        vec = vec[:, index]
        m[ID] = vec
        if ID in g1_id:
            m1[i1][:] = vec
            i1 += 1
        if ID in g2_id:
            m2[i2][:] = vec
            i2 += 1
        

    g1_mean = np.mean(m1, 0)
    g2_mean = np.mean(m2, 0)
    group1 = []
    group2 = []
    for ID, vec in read_vector:
        if ID not in m:
            continue
        s1 = np.sum( np.array(m[ID]) * g1_mean)
        s2 = np.sum( np.array(m[ID]) * g2_mean)
        if s1 > s2:
            group1.append( (ID, vec) )
        else:
            group2.append( (ID, vec) )
    return group1, group2 

def partition_reads(read_vector, g_id, MIN_GROUP_SIZE, ENTROPY_TH):
    print "partition group:", g_id
    if len(read_vector) < MIN_GROUP_SIZE:
        return (("%d" % g_id, read_vector), )
    
    smap = {"A":1,"C":1,"G":1,"T":1,"*":-1," ":0}
    
    m = np.zeros( (len(read_vector), len(read_vector[0][1])) )
    i = 0
    for ID, node in read_vector:
        m[i] = [smap[c] for c in node]
        i += 1

    pos_candidates = []
    ce = col_entropy(m)
    total_index = []
    for pos in range(m.shape[1]):
        entropy = ce[pos]
        total_index.append( (entropy, pos) )
        if entropy > ENTROPY_TH:
            pos_candidates.append( (entropy, pos) )
    pos_candidates = [x[1] for x in pos_candidates]

    total_index = range(m.shape[1])
    #total_index.sort()
    #total_index.reverse()
    #n_pos = min(64, len(pos_candidates) * 2)
    #total_index  = [x[1] for x in total_index[:n_pos]]
    #total_index.sort()
    print "number_of_candidate", len(total_index)

    partition_scores = []
    for pos in pos_candidates:
        partition_scores.append( ( (partition_score(m, pos, total_index), -ce[pos]), pos) )
    partition_scores.sort()
   
    if len(partition_scores) == 0:
        return (("%d" % g_id, read_vector, m), )
    print len(read_vector), partition_scores[0]
    print "------------------"
    if partition_scores[0][0][0] > 0:
        return (("%d" % g_id, read_vector, m), )
    
    #rtn = ["groups_%d" % g_id]
    rtn = []
    group1 = []
    group2 = []
    split_pos = partition_scores[0][1]
    for rid, vec in read_vector:
        if vec[split_pos] == "*":
            group1.append( (rid, vec) )
        elif vec[split_pos] in ("A", "C", "G", "T"):
            group2.append( (rid, vec) )
    group1, group2 = regroup(group1, group2, total_index)

    if len(group1) < MIN_GROUP_SIZE or len(group2) < MIN_GROUP_SIZE:
        return (("%d" % g_id, read_vector, m), )
    else:
        rtn.extend( partition_reads(group1, g_id * 2 + 1, MIN_GROUP_SIZE, ENTROPY_TH ) )
        rtn.extend( partition_reads(group2, g_id * 2 + 2, MIN_GROUP_SIZE, ENTROPY_TH ) )
    
    return rtn
    


def level2_partition(read_ids, read_file, ref_file, WORK_DIR, ENTROPY_TH, MIN_GROUP_SIZE):

    level2_group = []

    print "s:" + read_file
    ignore_indel = False
    tmp_cns = os.path.join( WORK_DIR, "tmp_cns.fa")
    aln_g, seq = get_consensus(read_file, ref_file, tmp_cns, "tmp_cns", 
                               ENTROPY_TH,
                               hp_correction = False, 
                               min_iteration = 4, 
                               max_num_reads = 150,
                               entropy_th = 0.65,
                               min_cov = 8,
                               max_cov = 200,
                               nproc = 16,
                               mark_lower_case = False,
                               use_read_id = False)
        
    aln_g = construct_aln_graph_from_fasta(read_file, 
                                            tmp_cns, 
                                            max_num_reads = 5000, 
                                            max_cov = 5000, 
                                            remove_in_del = True, 
                                            nproc=16, 
                                            use_read_id=True)

    seq, c_data = aln_g.generate_consensus(min_cov=0)


    if len(read_ids) < MIN_GROUP_SIZE:
        print "group element < %d, not splitting" % MIN_GROUP_SIZE
        level2_group.append( ( read_ids, seq, c_data, "-" ) )
        return level2_group


    rv, hen = read_node_vector(aln_g, ENTROPY_TH, entropy_th = ENTROPY_TH)
    if len(rv) == 0:
        print "not high entropy node cnd#1, not splitting"
        level2_group.append( ( read_ids, seq, c_data, "+" ) )
        return level2_group
    if len(rv.values()[0]) == 0:
        print "not high entropy node cnd#2, not splitting"
        level2_group.append( ( read_ids, seq, c_data, "+" ) )
        return level2_group

    aln_data = []
    for r in rv:
        aln_data.append( (r, "".join(rv[r])) )
    aln_data.sort(key=lambda x: x[1][0])
    aln_data.reverse()

    read_groups = partition_reads(aln_data, 0, MIN_GROUP_SIZE, ENTROPY_TH)
    print sum([len(rg[1]) for rg in read_groups]), [ ( rg[0], len(rg[1]) ) for rg in read_groups ]

    if len(read_groups) == 1:
        print "level 1 splitting fail, not splitting"
        level2_group.append( ( read_ids, seq, c_data, "+" ) )
        return level2_group

    group_id = 0
    for rg in read_groups:
        r_ids = set() 
        rd = rg[1]
        for r, d in rd:
            r_ids.add( r )

        f = FastaReader(read_file)
        tmp_hash = hash( "%s_%d" % (read_file, group_id) )
        read_out_file = os.path.join( WORK_DIR, "tmp_reads_%d.fa"  % tmp_hash )
        fetch_read(read_file, read_out_file, r_ids)

        tmp_cns = os.path.join( WORK_DIR, "tmp_cns_%s.fa" % tmp_hash )
        with open(tmp_cns,"w") as f:
            print >>f, ">tmp_cns"
            print >>f, seq

        group_id += 1
        level2_group.extend( level2_partition(r_ids, read_out_file, tmp_cns, WORK_DIR,
                                                     ENTROPY_TH, MIN_GROUP_SIZE, ) )
    return level2_group

class Clusense( object ):
    
    def __init__(self, read_file, ref_file, output_dir,
                                            threshold=None,
                                            entropy=None, 
                                            prefix=None, 
                                            min_group=None):
        self.read_file = read_file
        self.ref_file = ref_file
        self.output_dir = output_dir
        self.threshold = threshold
        self.entropy = entropy
        self.prefix = prefix
        self.min_group = min_group
        # Validate and run
        self.validate_args()
        self.initialize_logger()
        self.run()

    def validate_args(self):
        # Check the output directory, and create it if needed
        if self.output_dir is None:
            self.output_dir = args.read_file.split('.')[0] + '_wd'
        self.output_dir = os.path.abspath( self.output_dir )
        try:
            os.makedirs( self.output_dir )
        except:
            pass
        # Set the entropy value based on the supplied threshold
        if self.entropy is None:
            self.entropy = calculate_entropy( self.threshold )

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        print self.log

    def run(self):
        tmp_cns = os.path.join( self.output_dir, "tmp_cns.fa")
        aln_g, seq = get_consensus(self.read_file, 
                                   self.ref_file, 
                                   tmp_cns, 
                                   "tmp_cns",
                                   self.entropy,
                                   hp_correction = True, 
                                   min_iteration = 4, 
                                   max_num_reads = 150,
                                   entropy_th = 0.65,
                                   min_cov = 8,
                                   max_cov = 200,
                                   nproc = 16,
                                   mark_lower_case = False,
                                   use_read_id = False)

        aln_g = construct_aln_graph_from_fasta(self.read_file, 
                                               tmp_cns, 
                                               max_num_reads = 5000, 
                                               max_cov=5000, 
                                               remove_in_del = False, 
                                               nproc=16, 
                                               use_read_id = False)
        seq, c_data = aln_g.generate_consensus(min_cov=0, compute_qv_data= True)

        cns = os.path.join( self.output_dir, "group_root_cns.fa")
        with open(cns,"w") as f:
            print >>f, ">group_root_cns" 
            print >>f, seq
            
        score = os.path.join( self.output_dir, "group_root.score")
        with open(score,"w") as f:
            for i in range(len(seq)):
                print >>f, i, seq[i], " ".join([str(c) for c in c_data[i]]), 1.0*c_data[i][0]/(c_data[i][3]+1)

        r_ids = set()
        for r in FastaReader( self.read_file ):
            r_ids.add(r.name)
            
        level2_group = level2_partition(r_ids, self.read_file, self.ref_file, 
                                                               self.output_dir, 
                                                               self.entropy, 
                                                               self.min_group)
        s = 0
        print "-------------------"
        group_id = 1

        summary_f = open(os.path.join( self.output_dir, "summary.txt" ), "w")

        for id_set, seq, c_data, status in level2_group:
            out_read_file =  os.path.join( self.output_dir, "group_%02d.fa" % group_id )
            fetch_read(self.read_file, out_read_file, id_set)
            cns = os.path.join( self.output_dir, "group_%02d_cns.fa" % group_id )

            with open(cns,"w") as f:
                print >>f, ">group_%02d_cns" % group_id
                print >>f, seq

            aln_g = construct_aln_graph_from_fasta(out_read_file, 
                                                   cns, 
                                                   max_num_reads = 5000, 
                                                   max_cov=5000, 
                                                   remove_in_del = False, 
                                                   nproc=16, 
                                                   use_read_id = False)
            seq, c_data = aln_g.generate_consensus(min_cov=0, compute_qv_data= True)

            with open(cns,"w") as f:
                print >>f, ">%s_group_%02d_cns" % (self.prefix, group_id)
                print >>f, seq
                
            score = c_data

            score = os.path.join( self.output_dir, "group_%02d.score" % group_id )
            with open(score,"w") as f:
                for i in range(len(seq)):
                    print >>f, i, seq[i], " ".join([str(c) for c in c_data[i]]), 1.0*c_data[i][0]/(c_data[i][3]+1)

            print >> summary_f, "group_%02d" % group_id, len(id_set)
            s += len(id_set)
            group_id += 1
        print >> summary_f, "total", s
        summary_f.close()

if __name__ == "__main__":
    import argparse
    desc = "Graph-model based hapolytype read separation"
    parser = argparse.ArgumentParser(description=desc)

    add = parser.add_argument
    add("read_file", metavar="READS", 
        help="Fasta-format file of sequence reads to separate")
    add("reference", metavar="REFERENCE", 
        help="Fasta-format file of the reference sequence to use")
    add("-o", "--output_dir", default=None, 
        help="Name of the directory to output results to")
    add("-m", "--min_group", type=int, default=MIN_GROUP,
        help="Minimum group size to require from each cluster ({0})".format(MIN_GROUP))
    add("-p", "--prefix", default=PREFIX, 
        help="Prefix to prepend before the name of any consensus sequence ({0})".format(PREFIX))
    add("-t", "--threshold", type=float, default=THRESHOLD,
        help="Threhold value to use for calculating entropy cutoffs ({0})".format(THRESHOLD))
    add("-e", "--entropy", type=float, default=None,
        help=argparse.SUPPRESS)
    args = parser.parse_args()

    Clusense( args.read_file, args.reference, args.output_dir, 
              args.threshold, args.entropy, args.prefix, args.min_group)
