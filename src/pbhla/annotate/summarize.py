import os, logging

from pbcore.io.FastaIO import FastaReader

from pbhla.io.BlasrIO import BlasrReader
from pbhla.fasta.utils import fasta_size

log = logging.getLogger()

PCTID_THRESH = 97.0

class ContigInfo(object):
    def __init__(self, locus, contig, length, count, hit, pctid):
        self.locus = locus
        self.contig = contig
        self.length = int(length) if length.isdigit() else length
        self.count = int(count) if count.isdigit() else count
        self.hit = hit
        self.pctid = float(pctid)

    def __str__(self):
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(self.locus,
                                                     self.contig,
                                                     self.length,
                                                     self.count,
                                                     self.pctid,
                                                     self.hit)

def parse_contig_info( handle, locus ):
    try:
        info = [locus] + handle.next().strip().split()
    except:
        info = [locus, '-\t\t\t', '-', '-', '-', '0.0']
    return ContigInfo(*info)

def summarize_resequenced( locus_file, blasr_file, output_file ):
    result_dict = {}
    for entry in BlasrReader( blasr_file ):
        name = entry.qname.split('|')[0]
        name = name[:-4]
        result_dict[name] = [entry.pctsimilarity, entry.tname]
    with open( output_file, 'w') as output:
        for line in open( locus_file ):
            if line.startswith('Locus'):
                print >> output, 'Locus\tContig\tLength\tCount\tPctId\tHit\tReseqPctId\tReseqHit'
            else:
                parts = line.strip().split()
                name = parts[1]
                if name == '-': 
                    continue
                reseq_results = result_dict[name]
                output_parts = parts + reseq_results
                print >> output, '\t'.join( output_parts )

def meta_summarize_contigs( contig_files, output_file, excluded=[] ):
    with open(output_file, 'w') as output:
        print >> output, "Locus\tContig\tLength\tCount\tHit\tPctId"
        for filepath in sorted(contig_files):
            log.info( filepath )
            locus = filepath.split('_')[-2]
            log.info( locus )
            if locus in excluded:
                continue
            dummy = ContigInfo(locus, '-\t\t\t', '-', '-', '-', '0.0')
            with open(filepath, 'r') as handle:
                handle.next()
                first = parse_contig_info( handle, locus )
                second = parse_contig_info( handle, locus )
                while second.hit == first.hit:
                    second = parse_contig_info( handle, locus )
                if second.pctid < PCTID_THRESH:
                    second = dummy
                print >> output, first
                print >> output, second

def summarize_contigs( sequence_file, subread_fofn, locus_dict, blasr_file, output_dir):
    groups = group_by_locus( locus_dict )
    lengths = parse_fasta_lengths( sequence_file )
    sizes = parse_subread_counts( subread_fofn )
    hits = parse_blasr_alignment( blasr_file )
    # Summarize each locus group
    summary_outputs = []
    for locus, group in groups.iteritems():
        if locus == 'N/A':
            continue
        summary = summarize_group( group, lengths, sizes, hits )
        output = write_summary( locus, summary, output_dir )
        summary_outputs.append( output )
    return summary_outputs

def parse_blasr_alignment( blasr_file ):
    hits = {}
    for entry in BlasrReader( blasr_file ):
        name = entry.qname
        if name.endswith('_cns'):
            name = entry.qname[:-4]
        hits[name] = (entry.tname, entry.pctsimilarity)
    return hits

def parse_fasta_lengths( fasta_file ):
    lengths = {}
    for record in FastaReader( fasta_file ):
        name = record.name
        if name.endswith('_cns'):
            name = name[:-4]
        lengths[name] = len(record.sequence)
    return lengths

def parse_subread_counts( subread_fofn ):
    sizes = {}
    with open( subread_fofn ) as handle:
        for filepath in handle:
            filepath = filepath.strip()
            filename = os.path.basename( filepath )
            contig_name = filename.split('.')[0]
            if contig_name.startswith('Resequenced_'):
                contig_name = '_'.join( contig_name.split('_')[1:] )
            sizes[contig_name] = fasta_size( filepath )
    return sizes

def group_by_locus( locus_dict ):
    groups = {}
    for key, value in locus_dict.iteritems():
        try:
            groups[value].append( key )
        except:
            groups[value] = [ key ]
    return groups

def summarize_group( group, lengths, sizes, hits):
    summary = []
    for contig in group:
        if contig.endswith('_cns'):
            contig = contig[:-4]
        try:
            hit, pctid = hits[contig]
            summary.append( (contig, 
                             lengths[contig], 
                             sizes[contig],
                             hit,
                             pctid) )
        except:
            pass
    return sorted(summary, key=lambda x: x[2], reverse=True)

def write_summary( locus, summary, output_dir ):
    output_file = 'Locus_{0}_Summary.txt'.format( locus )
    output_path = os.path.join( output_dir, output_file )
    with open(output_path, 'w') as handle:
        handle.write( 'Contig\tLength\tCount\tBestHit\tPctId\n' )
        for item in summary:
            line = '\t'.join( map(str, item) ) 
            handle.write( line + '\n' )
    return output_path
