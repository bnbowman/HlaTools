import os, logging

from collections import namedtuple

from pbcore.io.FastaIO import FastaReader

from pbhla.io.BlasrIO import BlasrReader
from pbhla.fasta.utils import fasta_size

contig_info = namedtuple('contig_info', 'contig length reads hit pctid')

log = logging.getLogger()

def meta_summarize_contigs( contig_files, output_dir ):
    summaries = []
    contig_summary = os.path.join( output_dir, 'Sample_Calls.txt')
    with open(contig_summary, 'w') as output:
        for filepath in sorted(contig_files):
            locus = filepath.split('_')[-2]
            with open(filepath, 'r') as handle:
                handle.next()
                try:
                    first = [locus] + handle.next().strip().split()
                except:
                    first = [locus, '-', '-', '-', '-', '-']

                second = handle.next().strip().split()
                print [locus] + first
                print [locus] + second

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
