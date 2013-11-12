#! /usr/bin/env python

import sys
import re
from math import ceil

infile = sys.argv[1] ### the .txt file from IMGT
outfile = sys.argv[2] ### output file basename
### will output $basename.afa and $basename.info

### third argument is either 'gDNA' for parsing gDNA alignments or 'cDNA' for parsing cDNA alignments
### fourth argument is starting base index (determine this by doing `head` on the IMGT alignment. It is the first number you see.)
### firth argument is the starting codon index (for cDNA only, determine this by doing `head` on the cDNA IMGT alignment. It is the second number you see.)

class Alignment:
    def __init__(self, name, backbone, isref, pos_dict={}, features={}, codon_features={}, start_pos=0, start_codon=0):
	self.name = name
	self.currentpos=0
	self.canonpos=start_pos
	self.pos_dict=pos_dict
	self.sequence=[]
	self.backbone=backbone
	self.isref=isref
	self.features=features
	self.status=0
	self.codon_counter=-1
	self.current_codon=start_codon
	self.codon_features=codon_features
    def add(self, base):
	if base == "*":
	    self.sequence.append("-")
	    if self.isref:
		self.features[self.currentpos] = self.status
		self.backbone[self.currentpos] = base
		self.codon_features[self.currentpos]=self.current_codon
		self.pos_dict[self.currentpos] = self.canonpos
		self.canonpos+=1
	    self.currentpos+=1 
	elif base == ".":
	    self.sequence.append("-")
	    if self.isref:
		self.features[self.currentpos] = self.status
		self.backbone[self.currentpos] = base
		self.codon_features[self.currentpos]=self.current_codon
		self.pos_dict[self.currentpos] = self.canonpos
	    self.currentpos+=1; 
	elif base == "-":
	    if self.isref:
		self.sequence.append("-")
		self.features[self.currentpos] = self.status
		self.codon_features[self.currentpos]=self.current_codon
		self.backbone[self.currentpos] = base
		self.pos_dict[self.currentpos] = self.canonpos
		self.canonpos+=1
	    else:
		self.sequence.append(self.backbone[self.currentpos])
		
	    self.currentpos+=1
	elif base in ['A', 'G', 'C', 'T']:
	    self.sequence.append(base)
	    if self.isref:
		
		self.features[self.currentpos] = self.status
		self.codon_counter+=1
		if self.codon_counter == 3:
		    self.current_codon+=1	
		    if self.current_codon == 0:
			self.current_codon+=1
		    self.codon_counter=0
		self.codon_features[self.currentpos]=self.current_codon
		self.backbone[self.currentpos] = base
		self.pos_dict[self.currentpos] = self.canonpos
		self.canonpos+=1
	    self.currentpos+=1
	elif base == "|":
	    if self.isref:
		print self.currentpos
		print self.canonpos
	    self.status +=1	
	else:
	    pass
	if self.canonpos==0:
	    self.canonpos=1
	
    def addline(self, line):
	for character in line:
	    self.add(character)
    def write(self, f):
	print >>f, ">"+self.name
	print >>f, ''.join(self.sequence)
		
	
		

with open(infile, "r") as f:
    firstline = True
    backbone={}; alignments={}; pos_dict={}
    for line in f:
	if len(line.strip()) < 10:
	    continue
	line = re.split("  +", line.strip() )
	if len(line) != 2:
	    continue
	if len(line[1]) < 10:
	    continue
	allele = line[0].strip()
	line = line[1].strip()
	if firstline:
		isref = True	
	else:
		isref = False
	if allele not in alignments.keys():
	    alignments[allele] = Alignment(allele, backbone, isref, pos_dict = pos_dict, start_pos=int(sys.argv[4]), start_codon=int(sys.argv[5]) )
	    if isref:
		ref_alignment=alignments[allele]
	    alignments[allele].addline(line)
	else:
	    alignments[allele].addline(line)
	firstline=False

with open(outfile+'.afa', "w") as of:
    for alignment in alignments.itervalues():
	alignment.write(of)

feature_dict = ref_alignment.features
codon_dict=ref_alignment.codon_features
pos_dict=ref_alignment.pos_dict
feature_names={}
max_item=-1
for item in feature_dict.values():
    if sys.argv[3] == 'gDNA':
	if (item % 2) == 0:
	    feature_names[item]='intron'+str(int(item/2))
	elif (item % 2) != 0:
	    feature_names[item]='exon'+str(int(ceil(item/float(2))))
	if item > max_item: max_item = item
    else:
	feature_names[item]='exon'+str(int(item)+1)
if sys.argv[3] == 'gDNA':
    feature_names[0]='5UTR'
    feature_names[max_item]='3UTR'	
for item in feature_dict.keys():
    feature_dict[item] = feature_names[feature_dict[item]]
with open(outfile+'.info', "w") as of:
    for i in xrange(len(feature_dict.keys())):
	if sys.argv[3] == 'gDNA':
	    print>>of, '%s %s %s' % (i,pos_dict[i],feature_dict[i])
	elif sys.argv[3] == 'cDNA':
	    print>>of, '%s %s %s %s' % (i,pos_dict[i],feature_dict[i], codon_dict[i])	
	
