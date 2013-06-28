from pbcore.io.FastaIO import FastaReader
from collections import namedtuple

VALID_BASES = frozenset(['A', 'G', 'C', 'T'])

def MSA_aligner(query, target, name=None, mode = "string_to_fasta"):
    assert mode in ['string_to_fasta', 'string_to_string']
    
    if mode == 'string_to_fasta':
	assert name != None
	results=[]
	start, end = examine_profile(query)
	for r in FastaReader(target):
	    if r.name == name: 
	        continue
	    total = 0
	    match = 0
	    entry = {}
	    realpos = 1
	    for item in zip(query[start:end+1],r.sequence[start:end+1]):
		if item[0] == item[1]:
		    match+=1
		elif item[0] != item[1]:
		    entry[total+start]=(item[0], item[1], realpos)	
		if item[0] in VALID_BASES:
		    realpos+=1
		total+=1
	    results.append([ r.name,
	                     match,
	                     total,
	                     match/float(total), 
	                     entry, 
	                     query[start:end+1], 
	                     r.sequence[start:end+1]])
	results.sort(key = lambda x: x[3], reverse=True)
        return results[0][:5]

    if mode == 'string_to_string':
	compare_start1, compare_end1 = examine_profile(query)
	compare_start2, compare_end2 = examine_profile(target)
	compare_start = sorted([compare_start1, compare_start2])[-1]
	compare_end = sorted([compare_end1, compare_end2])[0]
	total=0; match=0; entry={}; qentry={}; tentry={}; qrealpos=1; trealpos=1
	for item in zip(query[compare_start:compare_end+1],target[compare_start:compare_end+1]):
	    if item[0] == item[1]:
		match+=1
	    elif item[0] != item[1]:
		entry[total+compare_start]=(item[0], item[1], qrealpos, trealpos)
		qentry[total+compare_start]=(item[0], item[1], qrealpos)
		tentry[total+compare_start]=(item[1], item[0], trealpos)
	    if item[0] in VALID_BASES:
		qrealpos+=1
	    if item[1] in VALID_BASES:
                trealpos+=1
	    total+=1
	info = namedtuple('info', 'qvars, tvars, vars, score, qaln, taln')
	results = info._make([ qentry, 
	                       tentry, 
	                       entry, 
	                       match/float(total), 
	                       query[compare_start:compare_end+1], 
	                       target[compare_start:compare_end+1] ])
	return results

def examine_profile(string):
    rstring = string[::-1]
    compare_start = 0
    compare_end = 0
    for base in rstring:
        if base in VALID_BASES:
            break	
        compare_end += 1
    compare_end = len(string)-(compare_end+1)
    for base in string:
        if base in VALID_BASES:
            break
        compare_start += 1
    return compare_start, compare_end
