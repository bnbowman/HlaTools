import os, logging

from pbcore.io.FastaIO import FastaReader
from pbhla.fasta.utils import fasta_size

log = logging.getLogger()

class ContigPicker( object ): 

    def __init__(self, sequence_file, subread_fofn, locus_dict, output_dir):
        self.sequence_file = sequence_file
        self.subread_fofn = subread_fofn
        self.locus_dict = locus_dict
        self.output_dir = output_dir
        self.run()

    def run(self):
        length_dict = parse_fasta_lengths( self.sequence_file )
        size_dict = parse_subread_counts( self.subread_fofn )
        locus_groups = group_by_locus( self.locus_dict )
        # Summarize each locus group
        for locus, group in locus_groups.iteritems():
            summary = summarize_group( group, length_dict, size_dict )
            write_summary( locus, summary, self.output_dir )


def parse_fasta_lengths( fasta_file ):
    lengths = {}
    for record in FastaReader( fasta_file ):
        lengths[record.name] = len(record.sequence)
    return lengths

def parse_subread_counts( subread_fofn ):
    sizes = {}
    with open( subread_fofn ) as handle:
        for filepath in handle:
            filepath = filepath.strip()
            filename = os.path.basename( filepath )
            contig_name = filename.split('.')[0]
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

def summarize_group( group, length_dict, size_dict):
    summary = []
    for contig in group:
        summary.append( (contig, length_dict[contig], size_dict[contig]) )
    return sorted(summary, key=lambda x: x[2], reverse=True)

def write_summary( locus, summary, output_dir ):
    output_file = 'Locus_{0}_Summary.txt'.format( locus )
    output_path = os.path.join( output_dir, output_file )
    with open(output_path, 'w') as handle:
        handle.write( 'Contig\tLength\tCount\n' )
        for item in summary:
            line = '\t'.join( map(str, item) ) 
            handle.write( line + '\n' )
