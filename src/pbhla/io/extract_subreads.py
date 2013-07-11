#! /usr/bin/env python

import os, random, logging

from pbcore.io.BasH5Reader import BasH5Reader
from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter

from pbhla.fasta.utils import write_fasta

MIN_LENGTH = 3000
MIN_SCORE = 0.8
MAX_COUNT = None

random.seed(42)
log = logging.getLogger(__name__)

def extract_subreads(input_file, output_file, min_length, min_score, max_count):
    """
    Extract, filter and subset subreads from Bas/Bax/Fofn Files
    """
    log.info('Extracting subreads from %s' % os.path.basename(input_file))
    log.debug('\tMinimum Length:\t%s' % min_length)
    log.debug('\tMinimum Score:\t%s' % min_score)
    log.debug('\tMax Count:\t%s' % max_count)

    subreads = []
    for i, filename in enumerate(_iterate_input_files( input_file )):
        if filename.endswith('.bas.h5') or filename.endswith('bax.h5'):
            subreads += _extract_from_bash5( filename, min_length, min_score )
        elif filename.endswith('.fa') or filename.endswith('.fasta'):
            subreads += _extract_from_fasta( filename, min_length )
    log.info("Extracted %s subreads from %s files" % (len(subreads), i+1))

    if max_count:
        subreads = _subset_subreads( subreads, max_count )

    write_fasta( subreads, output_file )

    log.info("Finished extracting subreads")

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

def _extract_from_bash5( bash5_file, min_length, min_score ):
    """
    Extract filtered subreads from a BasH5 or BaxH5 file
    """
    filename = os.path.basename( bash5_file )
    log.info("Extracting subreads from %s" % filename)

    records = []
    for zmw in BasH5Reader( bash5_file ):
        if zmw.readScore < min_score:
            continue
        for subread in zmw.subreads:
            if len(subread.basecalls()) < min_length:
                continue
            record = FastaRecord( subread.readName, subread.basecalls() )
            records.append( record )

    log.info('Found %s subreads that passed filters' % len(records))
    return records

def _extract_from_fasta( fasta_file, min_length ):
    """
    Extract filtered subreads from a Fasta file
    """
    filename = os.path.basename( fasta_file )
    log.info("Extracting subreads from %s" % filename)

    records = []
    for record in FastaReader( fasta_file ):
        if len(record.sequence) < min_length:
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
                      args.max_count )
