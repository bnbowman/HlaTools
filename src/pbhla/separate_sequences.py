import logging

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbhla.utils import get_base_sequence_name

log = logging.getLogger()

def separate_listed_sequences( fasta_file, good_values, good_output, bad_output ):
    """
    Separate a fasta file into two based on a supplied value list
    """
    with FastaWriter( good_output ) as good_handle:
        with FastaWriter( bad_output ) as bad_handle:
            for record in FastaReader( fasta_file ):
                name = get_base_sequence_name( record.name )
                if name in good_values:
                    good_handle.writeRecord( record )
                else:
                    bad_handle.writeRecord( record )

def separate_aligned_sequences( fasta_file, dictionary, good_values, good_output, bad_output ):
    """
    Separate a fasta file into two based on a supplied dictionary and value list
    """
    with FastaWriter( good_output ) as good_handle:
        with FastaWriter( bad_output ) as bad_handle:
            for record in FastaReader( fasta_file ):
                name = get_base_sequence_name( record.name )
                value = dictionary.get(name, "Unmapped")
                if value in good_values:
                    good_handle.writeRecord( record )
                else:
                    bad_handle.writeRecord( record )

def separate_sequences( fasta_file, dictionary, prefix='' ):
    """
    Separate a fasta file into multiple groups based on some dict
    """
    file_handles = {}
    for record in FastaReader( fasta_file ):
        name = get_base_sequence_name( record.name )
        group = dictionary.get( name, "Unmapped" )
        group_file = prefix + '_' + group + '.fasta'
        try:
            file_handles[group_file].writeRecord( record )
        except KeyError:
            file_handles[group_file] = FastaWriter( group_file )
            file_handles[group_file].writeRecord( record )
    return closed_file_handles( file_handles )

def closed_file_handles( file_handles ):
    """
    Does what is says on the tin
    """
    [file_handles[f].close() for f in file_handles]
    return [f for f in file_handles]