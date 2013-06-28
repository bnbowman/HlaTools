from pbcore.io.GffIO import Gff3Record

def create_annotation2(canonical_MSA_ref, new_MSA_ref, new_MSA_con, feature_map, seq_name, coord_dict=None):
    gff3_records=[]
    shift = 0
    realpos = 0
    attributes = {}
    for i in xrange(len(new_MSA_ref)):
	while True:
	    try:
		if new_MSA_ref[i] != canonical_MSA_ref[(i-shift)] and new_MSA_ref[i] == "-" and i != 0:
		    attributes['outside_canonical_coordinates'] = 1
		    shift += 1
		break
	    except IndexError:
		shift += 1
		continue
	if 'exon' in feature_map[(i-shift)].feature: 
	    isexon = 1 
	else: 
	    isexon = 0
	start = 0
	end = 0
	attributes['canonical_position'] = feature_map[(i-shift)].canonical_pos
	attributes['feature'] = feature_map[(i-shift)].feature
	attributes['codon'] = feature_map[(i-shift)].codon
	attributes['isexon'] = isexon
	attributes['basecall'] = new_MSA_con[i]
	if new_MSA_con[i] in ['A', 'G', 'C', 'T' ]:
	    realpos += 1
	    attributes['sequence_position'] = realpos
	    if coord_dict == None:
		start = realpos
		end = realpos
	    else:
		start = coord_dict[realpos]
		end = coord_dict[realpos]
	    record = Gff3Record(seq_name, start, end, 'base', attributes=attributes)
	    gff3_records.append(record)
    return gff3_records 

def create_annotation(canonical_MSA_ref, new_MSA_ref, new_MSA_con, coverage_map, feature_map, seq_name, quality_map, coord_dict=None):
    gff3_records=[]
    shift = 0
    realpos = 0
    for i in xrange(len(new_MSA_ref)):
	record = Gff3Record(seq_name)
	while True:
	    try:
		if new_MSA_ref[i] != canonical_MSA_ref[(i-shift)] and new_MSA_ref[i] == "-" and i != 0:
		    record.put('outside_canonical_coordinates', 1)
		    shift += 1
		break
	    except IndexError:
		shift += 1
		continue
	if 'exon' in feature_map[(i-shift)].feature: 
	    isexon = 1 
	else: 
	    isexon = 0
	record.start=0
	record.end=0
	record.put('canonical_position', feature_map[(i-shift)].canonical_pos )
	record.put('feature', feature_map[(i-shift)].feature )
	record.put('codon', feature_map[(i-shift)].codon )
	record.put('isexon', isexon )
	record.put('basecall', new_MSA_con[i] )
	if new_MSA_con[i] in ['A', 'G', 'C', 'T' ]:
	    realpos+=1
	    record.put('sequence_position', realpos )
	    if coord_dict == None:
		record.start=realpos
		record.end=realpos
	    else:
		record.start=coord_dict[realpos]
		record.end=coord_dict[realpos]
	    record.score=quality_map[realpos]
	    record.put('score', quality_map[realpos])
	    try:
		record.put('cov2', coverage_map[realpos] )
	    except:
		record.put('cov2', 'NA' )
	    record.type = 'base'
	    gff3_records.append(record)
    return gff3_records 

def create_var_annotation(consensus_seq, variant_map, name, coord_dict=None):
    ### output differences between sequence and best match into a variants.gff file
    lasttype=None; totalVar=''; gff3_records=[]
    for i in xrange(len(consensus_seq)):	
	try:
	    entry = variant_map[i]
	except KeyError:
	    lasttype = None
	    totalVar=''
	    continue
	reference = entry[1]
	variantSeq = entry[0]
	realpos = entry[2]
	assert variantSeq in ['A', 'G', 'C', 'T', '-' ]
	assert reference in ['A', 'G', 'C', 'T', '-' ]
	if variantSeq == "-" and reference in ['A', 'G', 'C', 'T' ]: varType = 'deletion'
	if reference == "-" and variantSeq in ['A', 'G', 'C', 'T' ]: varType = 'insertion'
	if reference in ['A', 'G', 'C', 'T' ] and variantSeq in ['A', 'G', 'C', 'T' ]: varType = 'substitution'	
	record = Gff3Record(name)
	record.type=varType 
	record.put('confidence', 50)
	record.put('coverage', 100)
	if varType == lasttype and lasttype == 'insertion':
	    totalVar+=variantSeq
	elif varType == lasttype and lasttype == 'deletion':
	    totalVar+=reference
	elif varType == 'insertion' and lasttype != varType:
	    totalVar=variantSeq
	elif varType == 'deletion' and lasttype != varType:
	    totalVar=reference
	if varType == 'insertion':
	    position = (realpos-len(totalVar))
	    if position < 1: position = 1
	    if coord_dict == None:
		new_pos = position
		record.end=position
		record.start=(position-len(totalVar))
	    else:
		new_pos = coord_dict[position]
		record.end=coord_dict[position]-1
                record.start=(coord_dict[position]-len(totalVar))
	    assert len(totalVar) >= 1
	    record.put('variantSeq', totalVar)
	    record.put('length', len(totalVar))
	    if varType == lasttype:
		gff3_records[-1] = record
	    else:
		gff3_records.append(record)
	elif varType == 'deletion':
	    position = realpos
	    if position < 1: position = 1
	    if coord_dict == None:
		new_pos = position
		record.start=position
		record.end=position
	    else:
		new_pos = coord_dict[position]
		record.start=coord_dict[position]
                record.end=coord_dict[position]
	    assert len(totalVar) >= 1
	    record.put('reference', totalVar)
	    record.put('length', len(totalVar))
	    if varType == lasttype:
		gff3_records[-1] = record
	    else:
		gff3_records.append(record)
	elif varType == 'substitution':
	    if coord_dict==None:
		record.start=(realpos)
		record.end=(realpos)
	    else:
		record.start=coord_dict[realpos]
                record.end=coord_dict[realpos]
	    record.put('variantSeq', variantSeq)
	    record.put('reference', reference)
	    record.put('length', 1)
	    gff3_records.append(record)	
	lasttype = varType	
    return gff3_records		
