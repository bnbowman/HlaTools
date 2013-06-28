from string import maketrans

from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter

from pbhla.utils import make_rand_string, getbash, runbash

COMPLEMENT = maketrans('ACGT', 'TGCA')

def reverse_complement( fasta_record ):
    rev_seq = fasta_record.sequence[::-1]
    rev_com_seq = rev_seq.translate( COMPLEMENT )
    return FastaRecord( fasta_record.name,
                        rev_com_seq )

def write_fasta(fasta_records, output_file):
    with FastaWriter( output_file ) as handle:
        for record in fasta_records:
            handle.writeRecord( record )

def fasta_size(fasta):
    try:
        f = FastaReader(fasta)
        count=0
        for read in f:
            count+=1
        return count
    except:
	return None

def fasta_length( fasta ):
    try:
        f = FastaReader( fasta )
    except:
        return 0
    return max([len(read.sequence) for read in f])

def extract_sequence(fasta, names):
    f = FastaReader(fasta)
    if isinstance(names, str):
        for r in f:
            if r.name == names:
                return r.sequence
    elif isinstance(names, list):
        output=[]
        for r in f:
            if r.name in names:
                output.append(r)
        return output

def copy_fasta( fasta, destination, name=None ):
    with FastaWriter( destination ) as handle:
        for record in FastaReader( fasta ):
            if name:
                record._name = name
            handle.writeRecord( record )

def combine_fasta( fasta_files, destination ):
    with FastaWriter( destination ) as handle:
        for fasta in fasta_files:
            for record in FastaReader( fasta ):
                handle.writeRecord( record )

"""
def trim_fasta(fasta_fn, outfile, mode="least_similar", criterion = "length", trim_to=2, pctid_threshold=float(0), coverage=None, aln_portion=float(0.5), max_chunk_size=int(20), tmpdir="/tmp", n_iter=10 ):
    tmpdir=tmpdir+"/"+make_rand_string()+"/"
    os.mkdir(tmpdir) 
    assert mode in ("least_similar", "redundant")
    assert criterion in ("length", "coverage")
    if mode == "least_similar":
        assert pctid_threshold == float(0)
    if criterion == "coverage":
        assert coverage != None
    f = FastaReader(fasta_fn)
    shuffled_reads=[]
    for r in f:
	shuffled_reads.append(r)
    random.shuffle(shuffled_reads)
    write_fasta(shuffled_reads, outfile)
    for iteration in xrange(n_iter):	
	read_number=0
	f = FastaReader(outfile)
	shuffled_reads=[]
	for r in f:
	    read_number+=1
	    shuffled_reads.append(r)
	random.shuffle(shuffled_reads)
	write_fasta(shuffled_reads, outfile)
	f = FastaReader(outfile)

	read_number=0; chunkn=0; chunk_files={}
	for r in f:
	    if read_number == max_chunk_size:
		read_number=0
		chunk_infile=chunk_outfile = outfile.split(".")[0]+"_inchunk%s.fasta" % (chunkn)
		chunk_infile=tmpdir+os.path.basename(chunk_infile)	
		chunk_outfile = outfile.split(".")[0]+"_chunk%s.fasta" % (chunkn)
		chunk_outfile=tmpdir+os.path.basename(chunk_outfile)
		runbash("mv %s %s; cp %s %s" % ( outfile, chunk_infile, chunk_infile, chunk_outfile ))
		chunk_files[chunk_infile]=chunk_outfile
		chunkn+=1
		read_number=0
	    if read_number==0:
		write_fasta([r], outfile, mode="w")
	    else:
		write_fasta([r], outfile, mode="a")
	    read_number+=1
	chunk_infile=chunk_outfile = outfile.split(".")[0]+"_inchunk%s.fasta" % (chunkn)
	chunk_infile=tmpdir+os.path.basename(chunk_infile)
	chunk_outfile = outfile.split(".")[0]+"_chunk%s.fasta" % (chunkn)
	chunk_outfile=tmpdir+os.path.basename(chunk_outfile)
	chunk_files[chunk_infile]=chunk_outfile
	runbash("mv %s %s; cp %s %s" % ( outfile, chunk_infile, chunk_infile, chunk_outfile ))

	for chunk_infile in chunk_files.keys():	
	    f = FastaReader(chunk_infile)
	    read_number=0
	    for r in f: read_number+=1
	    chunk_outfile = chunk_files[chunk_infile]
	    redundant_seqs=[]	
	    while 1:
		redundant_seq=None
		if mode == "least_similar" and read_number <= trim_to: break
		if read_number == 1: break
		pct_id_max=0
		nCandidates = read_number
		try:
		    self_alignment_output=parse_blasr(getbash("blasr %s %s -bestn %s -nCandidates %s -maxScore 0 -m 4 " \
			% (chunk_outfile, chunk_outfile, nCandidates, nCandidates) ), 4)
		except:
		    break
		for alignment in self_alignment_output:
		    if alignment.qname == alignment.tname: continue
		    t_aln_pct = (int(alignment.tend) - int(alignment.tstart))/float(alignment.tseqlength)
		    q_aln_pct = (int(alignment.qend) - int(alignment.qstart))/float(alignment.qseqlength)
		    ### is this an end to end alignment or just a partial one?
		    ### if its just partial skip it
		    full_aln = bool((t_aln_pct > aln_portion) or (q_aln_pct > aln_portion) )
		    if not full_aln: continue
		    if criterion == "length":
			if alignment.tseqlength > alignment.qseqlength:
			    to_remove = alignment.qname
			else:
			    to_remove = alignment.tname
		    if criterion == "coverage":
			if coverage[alignment.tname] > coverage[alignment.qname]:
			    to_remove = alignment.qname
			else:
			    to_remove = alignment.tname
		    if full_aln and float(alignment.pctid) > pctid_threshold:
			if float(alignment.pctid) > pct_id_max:
			    pct_id_max = float(alignment.pctid)
			    redundant_seq=to_remove

		if redundant_seq == None: break
		redundant_seqs.append(redundant_seq)

		f = FastaReader(chunk_infile)
		output_reads=[]
		for r in f:
		    if r.name not in redundant_seqs:
			output_reads.append(r)
		read_number = len(output_reads)
		write_fasta(output_reads, chunk_outfile, "w")
	    runbash("cat %s >> %s" % ( chunk_outfile, outfile) )
	for item in chunk_files.iteritems():
	    os.remove(item[0])
	    os.remove(item[1])
    os.rmdir(tmpdir)
    return 0
"""
