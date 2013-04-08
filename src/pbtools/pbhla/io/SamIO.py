#!/home/UNIXHOME/jquinn/HGAP_env/bin/python
from collections import namedtuple
 
class SamReader:
    def __init__(self, f):
        self.file = open(f, "r")
	self.entry = namedtuple('entry', 'qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, uid, score, aln_length')

    def __iter__(self):
        for line in self.file:
            if len(line)<1: continue
            if line[0]=='@': continue
	    line = line.strip().split()
	    score = float(line[13].split(":")[-1])
	    aln_length = int(float(line[16].split(":")[-1]))
            line = line[:12]
	    line.extend([score, aln_length])
            if len(line)==0: continue
	    yield self.entry._make(line)
