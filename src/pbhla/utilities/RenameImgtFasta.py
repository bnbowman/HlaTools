#! /usr/bin/env python

from pbcore.io import FastaReader, FastaWriter, FastaRecord

def rename_imgt_fasta( input_file, output_file ):
    with FastaWriter( output_file ) as handle:
        for record in FastaReader( input_file ):
            # Check that this is an IMGT-formatted FASTA record
            assert record.header.startswith('HLA:')

            # Extract the header and replace spaces with underscores
            new_header = record.header.strip().replace(' ', '_')

            # Create a new record with the same sequence and the type
            #    in place of it's id.
            new_record = FastaRecord(new_header, record.sequence)
            handle.writeRecord( new_record )

if __name__ == '__main__':
    import sys

    input_file = sys.argv[1]

    rename_imgt_fasta( input_file, sys.stdout )
