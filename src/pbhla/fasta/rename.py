#! /usr/bin/env python

import sys, logging

from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord
from pbhla.utils import read_dict_file, check_output_file

log = logging.getLogger(__name__)

def rename_fofn( input_fofn, output_fofn, name_key ):
    """
    Rename a FOFN of subread files
    """
    with open( output_fofn, 'w' ) as writer:
        with open( input_fofn, 'r' ) as handle:
            for line in handle:
                filename = line.strip()
                if filename:
                    renamed_file = filename.split('.')[0] + '_renamed.fasta'
                    rename_fasta( filename, renamed_file, name_key )
                    writer.write( renamed_file + '\n' )
    check_output_file( output_fofn )
    return output_fofn

def rename_fasta( input_file, output_file, name_key ):
    """
    Rename a single Fasta of subreads
    """
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
    check_output_file( output_file )
    return output_file

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    name_key = sys.argv[3]
    rename_fasta( input_file, output_file, name_key )
