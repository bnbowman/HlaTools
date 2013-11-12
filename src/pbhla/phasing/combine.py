#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

import os, logging

from pbcore.io.FastqIO import FastqReader, FastqWriter

log = logging.getLogger()

def combine_fastq( input_files, output_file):
    """
    Combine sequences from multiple Fastq files into one
    """
    log.info("Combining multiple Fastq outputs")
    record_counter = 0
    file_counter = 0
    with FastqWriter( output_file ) as writer:
        for filename in input_files:
            file_counter += 1
            for record in FastqReader( filename ):
                record_counter += 1
                writer.writeRecord( record )
    log.info("Found {0} consensus sequences in {1} outputs".format(record_counter,
                                                                   file_counter))
    return output_file

def combine_resequencing( input_dir, output_file ):
    """
    Combine all Resequencing results into a single Fastq file
    """
    log.info("Combining Resequencing outputs")
    record_counter = 0
    file_counter = 0
    with FastqWriter( output_file ) as writer:
        for result in find_resequencing_results(input_dir):
            file_counter += 1
            for record in FastqReader( result ):
                record_counter += 1
                writer.writeRecord( record )
    log.info("Found {0} consensus sequences in {1} outputs".format(record_counter,
                                                                   file_counter))
    return output_file

def find_resequencing_results( directory ):
    """
    Identify all Resequencing results in the supplied directory
    """
    for outer_entry in os.listdir( directory ):
        entry_path = os.path.join( directory, outer_entry )
        if os.path.isdir( entry_path ):
            for inner_entry in os.listdir( entry_path ):
                if inner_entry == 'consensus.fastq':
                    yield os.path.join( entry_path, inner_entry )

def combine_amp_analysis( input_dir, output_file ):
    """
    Combine all AmpAnalysis results into a single Fastq file
    """
    log.info("Combining AmpliconAnalysis outputs")
    record_counter = 0
    file_counter = 0
    with FastqWriter( output_file ) as writer:
        for result in find_amp_analysis_results(input_dir):
            file_counter += 1
            for record in FastqReader( result ):
                record_counter += 1
                writer.writeRecord( record )
    log.info("Found {0} consensus sequences in {1} outputs".format(record_counter,
                                                                   file_counter))
    return output_file

def find_amp_analysis_results( directory ):
    """
    Identify all AmpAnalysis results in the supplied directory
    """
    for outer_entry in os.listdir( directory ):
        entry_path = os.path.join( directory, outer_entry )
        if os.path.isdir( entry_path ):
            for inner_entry in os.listdir( entry_path ):
                if inner_entry == 'amplicon_analysis.fastq':
                    yield os.path.join( entry_path, inner_entry )