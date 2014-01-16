#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

from pbhla.io.BlasrIO import BlasrReader, BlasrWriter

def format_blasr_file( input_file, output_file ):
    with BlasrWriter( output_file ) as writer:
        with BlasrReader( input_file ) as reader:
            writer.write_header( reader.filetype )
            for record in reader:
                writer.write( record )

if __name__ == '__main__':
    import sys

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else sys.stdout

    format_blasr_file( input_file, output_file )