#! /usr/bin/env python

import os, logging

from pbhla.io.BlasrIO import BlasrReader
from pbhla.typing.HlaTypes import HlaType

log = logging.getLogger()

def summarize_typing( gDNA_align, cDNA_align, output=None ):
    """
    Pick the top N Amplicon Analysis consensus seqs from a Fasta by Nreads
    """
    if output is None:
        basename = '.'.join( gDNA_align.split('.')[:-1] )
        output = '%s.typing' % basename
    gDNA_data = _parse_alignment( gDNA_align )
    cDNA_data = _parse_alignment( cDNA_align )
    summaries = _summarize_hits( gDNA_data, cDNA_data )
    _write_type_summary( summaries, output )

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
        summary = _add_cDNA_summary( summary, cDNA_data[query] )
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
        if hasattr( hit, 'pctsimilarity' ):
            summary['gDNA_PctId'] = float( hit.pctsimilarity )
            summary['gDNA_Mismatch'] = 'N/A'
            summary['gDNA_Indel'] = 'N/A'
        elif hasattr( hit , 'nmat' ):
            similarity = float(hit.nmat) / (int(hit.nmat) + int(hit.nmis) + int(hit.nins) + int(hit.ndel))
            summary['gDNA_PctId'] = round(100*similarity, 2)
            summary['gDNA_Mismatch'] = int(hit.nmis)
            summary['gDNA_Indel'] = int(hit.nins) + int(hit.ndel)
        else:
            msg = 'Invalid alignment type!'
            log.error( msg )
            raise TypeError( msg )
    else:
        summary['Notes'].append('Multiple gDNA Matches')
        summary['gDNA_Type'] = 'N/A'
        summary['gDNA_PctId'] = 'N/A'
    return summary

def _add_cDNA_summary( summary, cDNA_hits ):
    """
    Summarize the cDNA hits for a consensus sequence
    """
    if len( cDNA_hits ) == 1:
        summary['cDNA_Type'] = HlaType.from_string( cDNA_hits[0].tname )
        summary['cDNA_PctId'] = float(cDNA_hits[0].pctsimilarity)
    elif len( cDNA_hits ) >= 2:
        summary['cDNA_Type'] = _combine_cDNA_hits( cDNA_hits )
        summary['cDNA_PctId'] = float(cDNA_hits[0].pctsimilarity)
    else:
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
    #elif ( summary['gDNA_Type'].gene   == summary['cDNA_Type'].gene   and
    #       summary['gDNA_Type'].field1 == summary['cDNA_Type'].field1 and
    #       summary['gDNA_Type'].field2 == summary['cDNA_Type'].field2 and
    #       summary['gDNA_Type'].field3 == summary['cDNA_Type'].field3):
    #    summary['Type'] = summary['gDNA_Type']
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

def _write_type_summary( summaries, output ):
    with open( output, 'w' ) as handle:
        handle.write('Sequence\tgDNA_Type\tgDNA_PctId\tgDNA_Mismatch\tgDNA_Indel\tcDNA_Type\tcDNA_PctId\tType\tNotes\n')
        for query, summary in summaries.iteritems():
            summary['Notes'] = ';'.join( summary['Notes'] )
            if summary['Notes']:
                summary['Notes'] = '"%s"' % summary['Notes']
            handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (summary['Sequence'],
                                                                   str(summary['gDNA_Type']).ljust(17),
                                                                   summary['gDNA_PctId'],
                                                                   summary['gDNA_Mismatch'],
                                                                   summary['gDNA_Indel'],
                                                                   str(summary['cDNA_Type']).ljust(14),
                                                                   summary['cDNA_PctId'],
                                                                   summary['Type'],
                                                                   summary['Notes']))

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    gDNA_align = sys.argv[1]
    cDNA_align = sys.argv[2]
    
    summarize_typing( gDNA_align, cDNA_align )
