#! /usr/bin/env python

import os, sys, random, logging

from pbcore.io.BasH5Reader import BasH5Reader
from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter
from pbhla.fasta.utils import write_fasta 
from pbhla.utils import check_output_file

MIN_LENGTH = 3000
MIN_SCORE = 0.8

log = logging.getLogger(__name__)

def extract_best_reads(input_file, 
                       output_file=None,
                       min_length=MIN_LENGTH, 
                       min_score=MIN_SCORE):
    """
    Extract, filter and subset subreads from Bas/Bax/Fofn Files
    """
    if output_file is None:
        basename = '.'.join( input_file.split('.')[:-1] )
        output_file = '%s.best.fasta' % basename
    log.info('Extracting subreads from %s' % os.path.basename(input_file))
    log.debug('\tMinimum Length:\t%s' % min_length)
    log.debug('\tMinimum Score:\t%s' % min_score)

    reads = []
    for i, filename in enumerate(_iterate_input_files( input_file )):
        reads += list( _extract_from_bash5( filename, min_length, min_score ))
    log.info("Extracted %s subreads from %s files" % (len(reads), i+1))

    write_fasta( reads, output_file )
    check_output_file( output_file )
    log.info("Finished extracting subreads")
    return output_file

def _iterate_input_files( input_file ):
    """
    Convert the input file to a list of absolute paths to Bash5s
    """
    if input_file.endswith('.bas.h5') or input_file.endswith('bax.h5'):
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
        zmwName = '%s/%s' % (zmw.baxH5.movieName, zmw.holeNumber)
        if zmw.readScore < min_score:
            continue
        #if zmw.ccsRead and len( zmw.ccsRead.basecalls() ) > min_length:
        #    yield FastaRecord( zmw.ccsRead.readName, zmw.ccsRead.basecalls() )
        #elif zmw.subreads:
        long_subreads = [s for s in zmw.subreads if len(s.basecalls()) > min_length]
        if len( long_subreads ) == 1:
            subread = long_subreads[0]
            yield FastaRecord( subread.readName, subread.basecalls() )
        elif len( long_subreads ) >= 2:
            ordered = sorted( long_subreads, key=lambda s: len(s.basecalls()), reverse=True )
            subread = ordered[0]
            yield FastaRecord( subread.readName, subread.basecalls() )
    log.info('Found %s subreads that passed filters' % len(records))

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
    args = parser.parse_args()

    extract_best_reads( args.input_file,
                        args.output,
                        args.min_length,
                        args.min_score )
