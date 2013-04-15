#!/home/UNIXHOME/jquinn/HGAP_env/bin/python
import re
import csv

from pbcore.io.FastaIO import FastaReader

class SubreadStats:
    def __init__(self, ref_fasta, locus_dict):
	self.total_reads=0
	self.total_fp_reads=0
	self.loci={}
	f = FastaReader(ref_fasta)
	for r in f:
	    ref_length = len(r.sequence)	
	    coverage = range(ref_length)
	    for i in range(ref_length): coverage[i]=0
	    self.loci[locus_dict[r.name]]={'aligned_reads' : 0,
		'aligned_fp_reads' : 0,
		'aligned_bp' : 0,
		'aligned_fp_bp' : 0,
		'raw_bp' : 0,
		'raw_fp_bp' : 0}

    @staticmethod
    def tabulate_subread_stats(data, name):
	avg_length = (data['raw_bp'])/float(data['aligned_reads'])
	avg_fp_length = (data['raw_fp_bp'])/float(data['aligned_fp_reads'])
	avg_aln_per_read = (data['aligned_bp'])/float(data['aligned_reads'])
	avg_aln_per_fp_read = (data['aligned_fp_bp'])/float(data['aligned_fp_reads'])
	return [name, data['aligned_reads'], data['aligned_fp_reads'], data['raw_bp'], data['raw_fp_bp'], \
		data['raw_bp'], avg_length, avg_fp_length, data['aligned_bp'], data['aligned_fp_bp'], avg_aln_per_read, avg_aln_per_fp_read]

    def write(self, outdir):
    	with open(outdir+"/subread_statistics.csv", "w") as of:
	    writer = csv.writer(of)
	    writer.writerow(["Locus","Reads","FP_Reads","BP","FP_BP","Avg_Length","Avg_FP_Length","Aln_BP","Aln_FP_BP","Avg_Aln_Per_Read","Avg_FP_Aln_Per_Read"])
	    all_stats = {'aligned_reads' : 0,
		    'aligned_fp_reads' : 0,
		    'aligned_bp' : 0,
		    'aligned_fp_bp' : 0,
		    'raw_bp' : 0,
		    'raw_fp_bp' : 0}
	    for item in self.loci.iteritems():
		    data=item[1]; name=item[0]
		    for item in data.iteritems(): 
			try:
			    all_stats[item[0]]+=item[1]
			except:
			    pass
		    writer.writerow(self.tabulate_subread_stats(data, name)) 
	    data = all_stats; name="All"
	    writer.writerow(self.tabulate_subread_stats(data, name))
