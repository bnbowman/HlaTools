import os, logging

from pbcore.io.FastaIO import FastaReader

from pbhla.io.BlasrIO import BlasrReader, BlasrM1, BlasrM5
from pbhla.fasta.utils import fasta_size
from pbhla.utils import get_base_sequence_name

log = logging.getLogger()

MIN_PCTID = 40.0
MIN_COUNT = 10

def summarize_typing( summary_file, gdna_file, cdna_file, output_file ):
    gdna_types = parse_typing( gdna_file )
    cdna_types = parse_typing( cdna_file )
    final_types = finalize_typing( gdna_types, cdna_types )
    combined_types = combine_typings( gdna_types, cdna_types, final_types )
    append_typing_results( summary_file, combined_types, output_file )

def append_typing_results( summary_file, combined_typings, output_file):
    typing_header = '\tGenType\tGenPctId\tExonType\tExonPctId\tType\t'
    with open(output_file, 'w') as output:
        with open(summary_file, 'r') as handle:
            header = handle.next().strip()
            output.write(header + typing_header)
            for line in handle:
                parts = line.strip().split()
                name = get_base_sequence_name( parts[1] )
                parts += combined_typings[name]
                output.write('\t'.join(parts) + '\n')

def combine_typings( gdna_types, cdna_types, final_types ):
    results = {}
    for name, final_type in final_types.iteritems():
        gdna_type, gdna_pctid = gdna_types[name]
        cdna_type, cdna_pctid = cdna_types[name]
        results[name] = [ gdna_type, 
                          gdna_pctid, 
                          cdna_type, 
                          cdna_pctid, 
                          final_type ]
    return results

def finalize_typing( gdna_types, cdna_types ):
    results = {}
    for name, gdna_typing in gdna_types.iteritems():
        gdna_type, gdna_pct = gdna_typing
        if name in cdna_types:
            cdna_type, cdna_pct = cdna_types[name]
            if gdna_type.startswith(cdna_type):
                results[name] = gdna_type
            else:
                results[name] = cdna_type
        else:
            results[name] = gdna_type
    return results

def parse_typing( typing_file ):
    results = {}
    with open( typing_file) as handle:
        for line in handle:
            if line.startswith('Locus'):
                continue
            parts = line.strip().split()
            name = get_base_sequence_name( parts[1] )
            typing = parts[2]
            pctid = parts[5]
            results[name] = [typing, pctid]
    return results

def _parse_blasr_alignment( blasr_file ):
    results = {}
    for entry in BlasrReader( blasr_file ):
        name = get_base_sequence_name( entry.qname )
        if isinstance(entry, BlasrM1):
            results[name] = [entry.tname, entry.pctsimilarity]
        elif isinstance(entry, BlasrM5):
            diffs = int(entry.nmis) + int(entry.nins) + int(entry.ndel)
            pctid = 100 * int(entry.nmat) / float(int(entry.nmat) + diffs)
            results[name] = [entry.tname, pctid]
    return results

def summarize_resequenced( locus_file, blasr_file, output_file ):
    blasr_typings = parse_blasr_alignment( blasr_file )
    with open( output_file, 'w') as output:
        for line in open( locus_file ):
            if line.startswith('Locus'):
                print >> output, 'Locus\tContig\tLength\tCount\tPctId\tHit\tReseqPctId\tReseqHit'
            else:
                parts = line.strip().split()
                name = parts[1]
                if name == '-': 
                    continue
                reseq_results = blasr_typings[name]
                output_parts = parts + reseq_results
                print >> output, '\t'.join( output_parts )

def summarize_contigs( sequence_file, subread_fofn, locus_dict, blasr_file, output_dir):
    groups = _group_by_locus( locus_dict )
    lengths = _parse_fasta_lengths( sequence_file )
    sizes = _parse_subread_counts( subread_fofn )
    hits = _parse_blasr_alignment( blasr_file )
    # Summarize each locus group
    summary_outputs = []
    for locus, group in groups.iteritems():
        if locus == 'N/A':
            continue
        fractions = _calculate_subread_fractions( group, sizes )
        summaries = _summarize_group( group, lengths, sizes, fractions, hits )
        filtered_summaries = list( _filter_summaries( summaries ))
        output = write_summary( locus, filtered_summaries, output_dir )
        summary_outputs.append( output )
    return summary_outputs

def _parse_fasta_lengths( fasta_file ):
    """
    Count the number of bases in each consensus sequence
    """
    lengths = {}
    for record in FastaReader( fasta_file ):
        name = record.name
        if name.endswith('_cns'):
            name = name[:-4]
        lengths[name] = len(record.sequence)
    return lengths

def _parse_subread_counts( subread_fofn ):
    """
    Count the number of subreads assocated with each consensus
    """
    sizes = {}
    with open( subread_fofn ) as handle:
        for filepath in handle:
            filepath = filepath.strip()
            filename = os.path.basename( filepath )
            contig_name = filename.split('.')[0]
            if contig_name.startswith('Allele_'):
                contig_name = '_'.join( contig_name.split('_')[1:] )
            if contig_name.startswith('Resequenced_'):
                contig_name = '_'.join( contig_name.split('_')[1:] )
            sizes[contig_name] = fasta_size( filepath )
    return sizes

def _calculate_subread_fractions( group, subread_counts ):
    """
    Calculate the fraction of subreads associated with each consensus
    """
    fractions = {}
    total = float(sum([subread_counts[key] for key in group]))
    for key, count in subread_counts.iteritems():
        fractions[key] = round(count / total, 3)
    return fractions

def _group_by_locus( locus_dict ):
    """
    Reverse the locus dictionary to group contigs by locus
    """
    groups = {}
    for key, value in locus_dict.iteritems():
        try:
            groups[value].append( key )
        except:
            groups[value] = [ key ]
    return groups

def _summarize_group( group, lengths, sizes, fractions, hits):
    """
    Iterate over the contigs in a group, summarizing each
    """
    summary = []
    for contig in group:
        if contig.endswith('_cns'):
            contig = contig[:-4]
        try:
            hit, pctid = hits[contig]
            summary.append( (contig, 
                             lengths[contig], 
                             sizes[contig],
                             fractions[contig],
                             hit,
                             pctid) )
        except:
            pass
    return sorted(summary, key=lambda x: x[2], reverse=True)

def _filter_summaries( summaries ):
    for summary in summaries:
        if summary[2] < MIN_COUNT:
            continue
        if float(summary[5]) < MIN_PCTID:
            continue
        yield summary

def write_summary( locus, summary, output_dir ):
    output_file = 'Locus_{0}_Summary.txt'.format( locus )
    output_path = os.path.join( output_dir, output_file )
    with open(output_path, 'w') as handle:
        handle.write( 'Contig\tLength\tCount\tFrac\tBestHit\tPctId\n' )
        for item in summary:
            line = '\t'.join( map(str, item) ) 
            handle.write( line + '\n' )
    return output_path
