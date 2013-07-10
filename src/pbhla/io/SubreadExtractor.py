#! /usr/bin/env python

import os, logging

from pbcore.io.BasH5Reader import BasH5Reader
from pbcore.io.FastaIO import FastaRecord, FastaWriter

MIN_LENGTH = 3000
MIN_SCORE = 0.80
DILUTION = 1.0

log = logging.getLogger()

def extract_subreads(input_file, output_file, min_length, min_score, dilution):
    log.info('Extracting subreads from "%s"' % input_file)
    log.debug('\tMinimum Length:\t%s' % min_length)
    log.debug('\tMinimum Score:\t%s' % min_score)
    log.debug('\tDilution Factor:\t%s' % dilution)

    file_list = _parse_input_filelist( input_file )
    log.info("Identified %s H5 files to extract subreads data from")
    writer = FastaWriter( output_file )

    for filename in file_list:
        count = 0
        for zmw in BasH5Reader( filename ):
            if zmw.readScore < min_score:
                continue
            for subread in zmw.subreads:
                if len(subread.basecalls()) < min_length:
                    continue
                count += 1
                record = FastaRecord( subread.readName, subread.basecalls() )
                writer.writeRecord( record )
        log.info('%s subreads passed filters from "%s"' % (count, filename))
    log.info("Finished extracting subreads\n")

"""
def extract_subreads(self, filename, fasta_handle, text_handle):
    for zmw in BasH5Reader( filename ):
        if zmw.readScore < self.min_score:
            continue
        adapter_starts = [z.readStart for z in zmw.adapters]
        adapter_ends = [z.readEnd for z in zmw.adapters]
        for subread in zmw.subreads:
            if len(subread.basecalls()) < self.min_length:
                continue
            record = FastaRecord( subread.readName, subread.basecalls() )
            fasta_handle.writeRecord( record )
            # If the read goes from adaptor to adaptor, record it as FP
            if ((subread.readStart in adapter_ends) and
                (subread.readEnd in adapter_starts)):
                text_handle.write( record.name + '\n' )
"""

def _parse_input_filelist( input_file ):
    """
    Convert the input file to a list of absolute paths to Bash5s
    """
    if input_file.endswith('.bas.h5') or input_file.endswith('bax.h5'):
        return [ input_file ]
    elif input_file.endswith('.fofn'):
        file_list = []
        with open( input_file, 'r' ) as handle:
            for line in handle:
                filepath = line.strip()
                file_list.append( filepath )
        return file_list
    else:
        msg = "Input file must be Bas.H5, Bax.H5 or FOFN!"
        log.info( msg )
        raise TypeError( msg )

if __name__ == '__main__':
    import argparse
    desc = 'Extract subreads and identify fullpass reads'
    parser = argparse.ArgumentParser( description=desc )

    add = parser.add_argument
    add('input_file', metavar='INPUT',
        help='Input Bas.H5 or FOFN file to extract reads from')
    add('output', metavar='PREFIX',
        help='Output prefix for both reads (Fasta) and Pass-Status (Txt)')
    add('--min_length', type=int, default=MIN_LENGTH,
        help='Minimum length for subreads ({0})'.format(MIN_LENGTH))
    add('--min_score', type=float, default=MIN_SCORE,
        help='Minimum ReadScore for subreads ({0})'.format(MIN_SCORE))
    add('--dilution', type=float, default=DILUTION,
        help='Fraction of reads to use ({0})'.format(DILUTION))
    args = parser.parse_args()

    extract_subreads( args.input_file,
                      args.output,
                      args.min_length,
                      args.min_score,
                      args.dilution )
