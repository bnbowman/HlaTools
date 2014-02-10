#! /usr/bin/env python

import os, sys, random, logging

from pbcore.io.BasH5IO import BasH5Reader
from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter

random.seed(42)
log = logging.getLogger(__name__)

def extract_subreads(input_file, 
                     output_file,
                     min_length,
                     max_length,
                     min_score,
                     min_snr,
                     max_count,
                     white_list=None):
    """
    Extract, filter and subset subreads from Bas/Bax/Fofn Files
    """
    log.info('Extracting subreads from %s' % os.path.basename(input_file))
    log.debug('\tMinimum Length:\t%s' % min_length)
    log.debug('\tMaximum Length:\t%s' % max_length)
    log.debug('\tMinimum Score:\t%s' % min_score)
    log.debug('\tMinimum SNR:\t%s' % min_snr)
    log.debug('\tMax Count:\t%s' % max_count)
    log.debug('\tWhitelisted ZMWs:\t%s' % white_list)

    if white_list:
        white_list = set( _parse_white_list( white_list ))

    subreads = []
    for i, filename in enumerate(_iterate_input_files( input_file )):
        if filename.endswith('.bas.h5') or filename.endswith('bax.h5'):
            subreads += _extract_from_bash5( filename, min_length, max_length, min_score, min_snr, white_list )
        elif filename.endswith('.fa') or filename.endswith('.fasta'):
            subreads += _extract_from_fasta( filename, min_length, max_length )
    log.info("Extracted %s subreads from %s files" % (len(subreads), i+1))

    if max_count:
        subreads = _subset_subreads( subreads, max_count )

    with FastaWriter( output_file ) as writer:
        for record in subreads:
            writer.writeRecord( record )

    log.info("Finished extracting subreads")
    return output_file

def _parse_white_list( white_list ):
    if white_list.endswith('.fasta') or white_list.endswith('.fa'):
        for record in FastaReader( white_list ):
            name = record.name.split()[0]
            zmw = '/'.join( name.split('/')[:2] )
            yield zmw
    elif white_list.endswith('.txt') or white_list.endswith('.ids'):
        with open( white_list ) as handle:
            for line in handle:
                name = line.strip().split()[0]
                zmw = '/'.join( name.split('/')[:2] )
                yield zmw

def _iterate_input_files( input_file ):
    """
    Convert the input file to a list of absolute paths to Bash5s
    """
    if input_file.endswith('.bas.h5') or input_file.endswith('bax.h5'):
        yield os.path.abspath( input_file )
    elif input_file.endswith('.fa') or input_file.endswith('.fasta'):
        yield os.path.abspath( input_file )
    elif input_file.endswith('.fofn'):
        file_list = []
        with open( input_file, 'r' ) as handle:
            for line in handle:
                filename = line.strip()
                yield os.path.abspath( filename )
    else:
        msg = "Input file must be Fasta, Bas.H5, Bax.H5 or FOFN!"
        log.info( msg )
        raise TypeError( msg )

def _extract_from_bash5( bash5_file, min_length, max_length, min_score, min_snr, white_list ):
    """
    Extract filtered subreads from a BasH5 or BaxH5 file
    """
    filename = os.path.basename( bash5_file )
    log.info("Extracting subreads from %s" % filename)

    records = []
    for zmw in BasH5Reader( bash5_file ):
        zmwName = '%s/%s' % (zmw.baxH5.movieName, zmw.holeNumber)
        if white_list and zmwName not in white_list:
            continue
        if zmw.readScore < min_score:
            continue
        if min(zmw.zmwMetric('HQRegionSNR')) < min_snr:
            continue
        for subread in zmw.subreads:
            if len(subread) < min_length:
                continue
            if len(subread) > max_length:
                continue
            record = FastaRecord( subread.readName, subread.basecalls() )
            records.append( record )

    log.info('Found %s subreads that passed filters' % len(records))
    return records

def _extract_from_fasta( fasta_file, min_length, max_length ):
    """
    Extract filtered subreads from a Fasta file
    """
    filename = os.path.basename( fasta_file )
    log.info("Extracting subreads from %s" % filename)

    records = []
    for record in FastaReader( fasta_file ):
        if len(record.sequence) < min_length:
            continue
        if len(record.sequence) > max_length:
            continue
        records.append( record )

    log.info("Found %s subreads that passed filters" % len(records))
    return records

def _subset_subreads( subreads, count ):
    """
    Select a subset of the subreads to return
    """
    if len(subreads) <= count:
        subset = subreads
    else:
        subset = random.sample( subreads, count )
    log.info('Selected %s subreads of %s' % (len(subset), len(subreads)))
    return subset

if __name__ == '__main__':
    import argparse
    desc = 'Extract subreads and identify fullpass reads'
    parser = argparse.ArgumentParser( description=desc )

    add = parser.add_argument
    add('input_file', 
        metavar='FILE',
        help='Input Bas.H5 or FOFN file to extract reads from')
    add('--output', 
        metavar='FILE', 
        default=sys.stdout,
        help='Output prefix for both reads (Fasta) and Pass-Status (Txt)')
    add('--min_length',
        metavar='INT',
        type=int, 
        default=MIN_LENGTH,
        help='Minimum length for subreads ({0})'.format(MIN_LENGTH))
    add('--min_score',
        metavar='FLOAT',
        type=float, 
        default=MIN_SCORE,
        help='Minimum ReadScore for subreads ({0})'.format(MIN_SCORE))
    add('--white_list',
        help='White list of ZMWs to extract reads from')
    add('--max_count', 
        metavar='INT',
        type=int, 
        default=MAX_COUNT,
        help='Maximum number of subreads to return ({0})'.format(MAX_COUNT))
    args = parser.parse_args()

    extract_subreads( args.input_file,
                      args.output,
                      args.min_length,
                      args.min_score,
                      args.white_list,
                      args.max_count )
