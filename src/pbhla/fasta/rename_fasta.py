#! /usr/bin/env python

import sys, logging

from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord
from pbhla.utils import read_dict_file

log = logging.getLogger(__name__)

def rename_fasta( input_file, output_file, name_key ):
    renaming_dict = read_dict_file( name_key )
    with FastaWriter( output_file ) as writer:
        for record in FastaReader( input_file ):
            old_name = record.name.split()[0]
            try:
                new_name = renaming_dict[old_name]
            except KeyError:
                msg = "Sequence name not found!"
                log.error( msg )
                raise KeyError( msg )
            new_record = FastaRecord( new_name, record.sequence )
            writer.writeRecord( new_record )

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    name_key = sys.argv[3]
    rename_fasta( input_file, output_file, name_key )
