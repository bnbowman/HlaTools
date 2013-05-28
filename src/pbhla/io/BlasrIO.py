from collections import namedtuple

def parse_blasr(output, mode, strip_query_names = True):
    parsed_output=[]
    if mode == 1:
        entry = namedtuple('entry', 'qname, tname, qstrand, tstrand, score, pctsimilarity, tstart, tend, tlength, qstart, qend, qlength, ncells')
    if mode == 4:
        entry = namedtuple('entry', 'qname, tname, score, pctsimilarity, qstrand, qstart, qend, qseqlength, tstrand, tstart, tend, tseqlength, mapqv, ncells, clusterScore, probscore, numSigClusters')
    if mode == 5:
        entry = namedtuple('entry', 'qname, qlength, z1, qalength, qstrand, tname, tlength, z2, talength, tstrand, score, nmis, nins, ndel, zscore, qseq, matchvector, tseq')

    output = output.strip().split("\n")
    output = [ x.split() for x in output ]
    if strip_query_names:
	new_output=[]
	for line in output:
	    line[0] = line[0].split("/")[0]
	    line[0] = line[0].split("|")[0]
	    line[1] = line[1].split("|")[0]
	    new_output.append(line)
	output = new_output

    for line in output:
	alignment = entry._make(line)
        parsed_output.append(alignment)
    return parsed_output
