#! /usr/bin/env python

import os, logging
from operator import itemgetter

from pbhla.io.BlasrIO import BlasrReader
from pbhla.typing.HlaTypes import HlaType
from pbhla.utils import check_output_file

log = logging.getLogger()

def summarize_typing( gDNA_align, cDNA_align, output=None ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    if output is None:
        basename = '.'.join( gDNA_align.split('.')[:-1] )
        output = '%s.typing' % basename
    log.info("Writing summary of typing information to %s" % output)
    gDNA_data = _parse_alignment( gDNA_align )
    cDNA_data = _parse_alignment( cDNA_align )
    summaries = _summarize_hits( gDNA_data, cDNA_data )
    _write_type_summary( summaries, output )
    check_output_file( output )
    return output

def _parse_alignment( alignment ):
    """
    Parse the genomic typeings from the gDNA alignment
    """
    hits = {}
    for record in BlasrReader( alignment ):
        try:
            hits[record.qname].append( record )
        except:
            hits[record.qname] = [ record ]
    return hits

def _summarize_hits( gDNA_data, cDNA_data ):
    """
    """
    summaries = {}
    for query in gDNA_data:
        summary = { 'Sequence': query,
                    'Notes': []}
        summary = _add_gDNA_summary( summary, gDNA_data[query] )
        summary = _add_cDNA_summary( summary, cDNA_data, query )
        summary = _add_consensus_type( summary )
        summaries[query] = summary
    return summaries

def _add_gDNA_summary( summary, gDNA_hits ):
    """
    Summarize the gDNA hits for a consensus sequence
    """
    if len( gDNA_hits ) == 1:
        hit = gDNA_hits[0]
        summary['gDNA_Type'] = HlaType.from_string( hit.tname )
        summary['gLen'] = hit.qlength
        if hasattr( hit, 'pctsimilarity' ):
            summary['gDNA_PctId'] = round(float( hit.pctsimilarity ), 2)
            summary['gDNA_Mismatch'] = 'N/A'
            summary['gDNA_Indel'] = 'N/A'
        elif hasattr( hit , 'nmat' ):
            similarity = float(hit.nmat) / (int(hit.nmat) + int(hit.nmis) + int(hit.nins) + int(hit.ndel))
            summary['gDNA_PctId'] = round(100*similarity, 2)
            summary['gDNA_Mismatch'] = str(hit.nmis)
            summary['gDNA_Indel'] = str(int(hit.nins) + int(hit.ndel))
        else:
            msg = 'Invalid alignment type!'
            log.error( msg )
            raise TypeError( msg )
    else:
        summary['Notes'].append('Multiple gDNA Matches')
        summary['gDNA_Type'] = 'N/A'
        summary['gDNA_PctId'] = 'N/A'
    return summary

def _add_cDNA_summary( summary, cDNA_data, query ):
    """
    Summarize the cDNA hits for a consensus sequence
    """
    if query in cDNA_data:
        cDNA_hits = cDNA_data[query]
        summary['cLen'] = cDNA_hits[0].qlength
        if len( cDNA_hits ) == 1:
            summary['cDNA_Type'] = HlaType.from_string( cDNA_hits[0].tname ).cDNA_type
            summary['cDNA_PctId'] = round(float(cDNA_hits[0].pctsimilarity), 2)
        elif len( cDNA_hits ) >= 2:
            summary['cDNA_Type'] = _combine_cDNA_hits( cDNA_hits )
            summary['cDNA_PctId'] = round(float(cDNA_hits[0].pctsimilarity), 2)
        else:
            summary['cDNA_Type'] = 'N/A'
            summary['cDNA_PctId'] = 'N/A'
    else:
        summary['cLen'] = 'N/A'
        summary['cDNA_Type'] = 'N/A'
        summary['cDNA_PctId'] = 'N/A'
    return summary

def _add_consensus_type( summary ):
    """
    Combine the gDNA and cDNA results into a consensus type
    """
    if summary['gDNA_Type'] == summary['cDNA_Type']:
        summary['Type'] = summary['gDNA_Type']
    elif summary['gDNA_PctId'] == 100.0:
        summary['Type'] = summary['gDNA_Type']
    elif summary['gDNA_PctId'] != 100.0 and summary['cDNA_PctId'] == 100.0:
        summary['Type'] = summary['cDNA_Type']
        summary['Notes'].append( 'Possible novel genomic sequence' )
    else:
        summary['Type'] = 'N/A'
    return summary

def _combine_cDNA_hits( cDNA_hits ):
    """
    Combine a series of cDNA hits at the 3-digit level
    """
    types = [HlaType.from_string( hit.tname ) for hit in cDNA_hits]
    genes = list( set( [t.gene for t in types] ))
    field1s = list( set( [t.field1 for t in types] ))
    field2s = list( set( [t.field2 for t in types] ))
    field3s = list( set( [t.field3 for t in types] ))
    if len(field3s) == 1:
        return HlaType( gene=genes[0], field1=field1s[0],
                                       field2=field2s[0],
                                       field3=field3s[0] )
    elif len(field2s) == 1:
        return HlaType( gene=genes[0], field1=field1s[0], 
                                       field2=field2s[0] )
    elif len(field1s) == 1:
        return HlaType( gene=genes[0], field1=field1s[0] )
    return 'N/A'

def _write_header( handle, seq_len=45, gtype_len=21, ctype_len=18 ):
    handle.write('Sequence'.ljust(seq_len))
    handle.write('gLen'.ljust(6))
    handle.write('gType'.ljust(gtype_len))
    handle.write('gPctId'.ljust(7))
    handle.write('nMis'.ljust(5))
    handle.write('Indel'.ljust(6))
    handle.write('cLen'.ljust(5))
    handle.write('cType'.ljust(ctype_len))
    handle.write('cPctId'.ljust(7))
    handle.write('Type')
    handle.write('\n')

def _write_summary_line( handle, summary, seq_len=45, gtype_len=21, ctype_len=18 ):
    handle.write(summary['Sequence'].ljust(seq_len))
    handle.write(summary['gLen'].ljust(6))
    handle.write(str(summary['gDNA_Type']).ljust(gtype_len))
    handle.write(str(summary['gDNA_PctId']).ljust(7))
    handle.write(summary['gDNA_Mismatch'].ljust(5))
    handle.write(summary['gDNA_Indel'].ljust(6))
    handle.write(summary['cLen'].ljust(5))
    handle.write(str(summary['cDNA_Type']).ljust(ctype_len))
    handle.write(str(summary['cDNA_PctId']).ljust(7))
    handle.write(str(summary['Type']))
    handle.write('\n')

def _get_max_field_size( summaries, field_name ):
    fields = [s[field_name] for s in summaries.itervalues()]
    return max([len(str(f)) for f in fields]) + 1

def _write_type_summary( summaries, output ):

    # First we iterate over the various fields to find their max sizes
    seq_len = _get_max_field_size( summaries, 'Sequence')
    gtype_len = _get_max_field_size( summaries, 'gDNA_Type')
    ctype_len = _get_max_field_size( summaries, 'cDNA_Type')

    # Using those lengths, we write out the records to file
    with open( output, 'w' ) as handle:
        _write_header( handle, seq_len, gtype_len, ctype_len )
        keys = {k: str(v['gDNA_Type']) for k, v in summaries.iteritems()}
        for key, g_type in sorted(keys.iteritems(), key=itemgetter(1)):
            summary = summaries[key]
            _write_summary_line( handle, summary, seq_len, gtype_len, ctype_len )

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    gDNA_align = sys.argv[1]
    cDNA_align = sys.argv[2]
    
    summarize_typing( gDNA_align, cDNA_align )
