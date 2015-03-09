#!/home/UNIXHOME/bbowman/bin/python2.7

import sys
import re

from Bio import SeqIO

for filename in sys.argv[1:]:
    assert re.search('(\.fa|\.fsa|\.fasta)$', filename)
    for record in SeqIO.parse(filename, 'fasta'):
        print ">%s" % record.id.strip()
        print str(record.seq)
