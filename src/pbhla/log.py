__author__ = 'bbowman@pacificbiosciences.com'

import sys
import logging

LOG_FORMAT = "%(asctime)s [%(levelname)s - %(module)s] %(message)s"
TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
FORMATTER = logging.Formatter( LOG_FORMAT, TIME_FORMAT )

def add_stream_handler( logger, stream=sys.stdout, log_level=logging.INFO ):
    # Set up a simple Stream handler
    stream_handler = logging.StreamHandler( stream=stream )
    stream_handler.setFormatter( FORMATTER )
    stream_handler.setLevel( log_level )
    logger.addHandler( stream_handler )

def add_file_handler( logger, log_file='hla_pipeline.log', log_level=logging.INFO ):
    # Set a second handler for the log file
    file_handler = logging.FileHandler( log_file )
    file_handler.setFormatter( FORMATTER )
    file_handler.setLevel( log_level )
    logger.addHandler( file_handler )

def initialize_logger( logger, stream=None, log_file=None, debug=False):
    if debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logger.setLevel( log_level )

    if stream:
        add_stream_handler( logger, stream=stream, log_level=log_level )
    else:
        add_stream_handler( logger, log_level=log_level )

    if log_file:
        add_file_handler( logger, log_file=log_file, log_level=log_level )
    else:
        add_file_handler( logger, log_level=log_level )
    return logger