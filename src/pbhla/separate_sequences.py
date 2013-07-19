import os, logging

from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbhla.utils import get_base_sequence_name

log = logging.getLogger()

def parse_id_list( list_file ):
    log.info('Parsing list of sequence Ids')
    sequence_list = []
    with open(list_file) as handle:
        for line in handle:
            sequence_list.append( line.strip() )
    return set(sequence_list)

def parse_value_list( value_list ):
    if isinstance(value_list, str):
        return parse_id_list( value_list )
    elif isinstnace(value_list, list):
        return set(value_list)
    else:
        raise TypeError

def separate_listed_sequences( fasta_file, value_list ):
    log.info('Separating sequences by Id')
    id_list = parse_value_list( value_list )
    # Separate sequences is the sequence list
    sequence_dict = {}
    for record in FastaReader( fasta_file ):
        name = record.name.split()[0]
        if name in id_list:
            add_record( sequence_dict, record, 'selected' )
        else:
            add_record( sequence_dict, record, 'not_selected' )
    return sequence_dict

def separate_aligned_sequences( fasta_file, dictionary, value_list ):
    log.info('Separating sequences by alignment hit')
    id_list = parse_value_list( value_list )
    # Separate sequences is the sequence list
    sequence_dict = {}
    for record in FastaReader( fasta_file ):
        name = record.name.split()[0]
        value = get_value(dictionary, name)
        if value in id_list:
            add_record( sequence_dict, record, 'selected' )
        else:
            add_record( sequence_dict, record, 'not_selected' )
    return sequence_dict

def separate_sequences( fasta_file, dictionary=None ):
    log.info('Separating out individual sequences')
    sequence_dict = {}
    for record in FastaReader( fasta_file ):
        name = get_base_sequence_name( record.name )
        if dictionary:
            value = get_value( dictionary, name )
            add_record( sequence_dict, record, value )
        else:
            add_record( sequence_dict, record, name )
    return sequence_dict

def add_record( dictionary, record, key):
    try:
        dictionary[key].append( record )
    except:
        dictionary[key] = [ record ]

def get_value(dictionary, key):
    try:
        value = dictionary[key]
    except:
        value = 'Unmapped'
    if value == 'N/A':
        value = 'Unmapped'
    return value

def write_group(dictionary, key, output_file):
    log.info('Writing "{0}" sequences to "{1}"'.format(key, os.path.basename(output_file)))
    with FastaWriter( output_file ) as handle:
        try:
            for record in dictionary[key]:
                handle.writeRecord( record )
        except:
            log.warn('No records found associated with "%s"' % key)

def write_all_groups(dictionary, prefix):
    log.info('Writing fasta sequences from all loci to file')
    output_files = []
    for key in sorted(dictionary):
        output_file = '{0}_{1}.fasta'.format(prefix, key)
        write_group(dictionary, key, output_file)
        output_files.append( output_file )
    log.info('Finished writing sequences from all loci')
    return output_files
