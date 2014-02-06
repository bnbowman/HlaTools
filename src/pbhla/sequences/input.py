#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os, re
import logging, logging.config

from pbcore.io import FastqReader, FastqWriter

from pbhla import __LOG__
from pbhla.utils import valid_file, is_fasta, is_fastq

logging.config.fileConfig( __LOG__ )
log = logging.getLogger(__name__)

def get_input_file( input ):
    if os.path.isdir( input ):
        log.info("Input appears to be a directory, looking for sequence files")
        return get_amplicon_analysis_output( input )
    if valid_file( input ):
        if is_fasta( input ):
            log.info("Input appears to be a valid Fasta file")
            return input
        elif is_fastq( input ):
            log.info("Input appears to be a valid Fastq file")
            return input
        else:
            msg = "Input is not a valid Fasta or Fastq file!"
            log.error( msg )
            raise IOError( msg )
    else:
        msg = "Supplied input does not appear to be a valid AmpliconAnalysis file or directory"
        log.error( msg )
        raise IOError( msg )

def get_amplicon_analysis_output( directory ):
    if 'amplicon_analysis.fastq' not in os.listdir( directory ):
        msg = "Directory does not appear to contain Amplicon analysis output!"
        log.error( msg )
        raise IOError( msg )
    elif 'amplicon_analysis_chimeras_noise.fastq' in os.listdir( directory ):
        msg = 'Chimera/Noise sequence file detected, combining with primary output'
        log.info( msg )
        return combine_amplicon_analysis_files( directory )
    else:
        msg = 'No Chimera/Noise sequence file detected, using primary output file'
        log.info( msg )
        return os.path.join( directory, 'amplicon_analysis.fastq' )

def combine_amplicon_analysis_files( directory ):
    output_file = os.path.join( directory, 'amplicon_analysis.all.fastq' )
    with FastqWriter( output_file ) as handle:
        for input_file in ['amplicon_analysis.fastq', 'amplicon_analysis_chimeras_noise.fastq']:
            input_path = os.path.join( directory, input_file )
            for record in FastqReader( input_path ):
                handle.writeRecord( record )
    return output_file