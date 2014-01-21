__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
import re
import logging

PB_REGEX = 'm\d{6}_\d{6}_[a-zA-Z0-9]{4,6}_c\d{33}_s\d_p\d'

log = logging.getLogger(__name__)

def invalid_fasta_names( fasta_file ):
    for record in FastaReader( fasta_file ):
        name = record.name.split()[0]
        if not re.match(PB_REGEX, name):
            return True
    return False

def invalid_dict_names( name_dict ):
    if name_dict is None:
        return True
    log.info('Checking name dictionary for valid well names')
    for key, value in name_dict.iteritems():
        if not re.match(PB_REGEX, value):
            return True
    return False

def is_exe( file_path ):
    if file_path is None:
        return False
    return os.path.isfile(file_path) and os.access(file_path, os.X_OK)

def which(program):
    """
    Find and return path to local executables
    """
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def create_directory( directory ):
    # Skip if the directory exists
    if os.path.isdir( directory ):
        return
    try: # Otherwise attempt to create it
        os.mkdir( directory )
    except:
        msg = 'Could not create directory "{0}"'.format(directory)
        log.info( msg )
        raise IOError( msg )

def read_dict_file( dict_file ):
    dict_contents = {}
    with open(dict_file, 'r') as handle:
        for line in handle:
            try:
                key, value = line.strip().split()
                dict_contents[key] = value
            except:
                pass
    return dict_contents
