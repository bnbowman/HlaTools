import os, logging

from pbhla.io.BlasrIO import BlasrReader

log = logging.getLogger()

def is_fasta( filename ):
    if filename.endswith('.fasta') or filename.endswith('.fa'):
        return True
    return False

def is_fastq( filename ):
    if filename.endswith('.fastq') or filename.endswith('.fq'):
        return True
    return False

def count_hits( filename ):
    return len( list( BlasrReader( filename )))

def get_base_sequence_name( name ):
    name = name.split()[0]
    if name.endswith('|quiver'):
        name = name.split('|')[0]
    if name.endswith('_cns'):
        name = name[:-4]
    return name

def get_file_type( filename ):
    if filename.endswith('.fa') or filename.endswith('.fasta'):
        return 'fasta'
    elif filename.endswith('.fq') or filename.endswith('.fastq'):
        return 'fastq'
    elif filename.endswith('.fofn'):
        return 'fofn'
    elif filename.endswith('.bas.h5') or filename.endswith('.bax.h5'):
        return 'bas.h5'
    else:
        msg = 'File is not of a recognized filetype'
        log.error( msg )
        raise TypeError( msg )

def memoize(function):
    cache = {}
    def decorated_function(*args):
        if args in cache:
            return cache[args]
        else:
            val = function(*args)
            cache[args] = val
            return val
    return decorated_function

def cleanup_directory( directory ):
    for entry in os.listdir( directory ):
        removal_flag = False
        if entry.endswith('aln') or entry.endswith('aln_unsorted'):
            removal_flag = True
        if entry.startswith('tmp_cns_') or entry.startswith('tmp_reads_'):
            removal_flag = True
        if removal_flag:
            try:
                os.remove( os.path.join( directory, entry) )
            except:
                pass

def write_list_file( file_list, output_file ):
    with open(output_file, 'w') as handle:
        for filename in file_list:
            print >> handle, filename

def read_list_file( list_file ):
    list_contents = []
    with open(list_file, 'r') as handle:
        for line in handle:
            value = line.strip().split()[0]
            if value:
                list_contents.append( value )
    return list_contents

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

def cross_ref_dict( query_dict, ref_dict ):
    new_dict = {}
    for key in query_dict:
        old_value = query_dict[key]
        if old_value.startswith('HLA:'):
            old_value = old_value.split('_')[0]
        try:
            new_value = ref_dict[old_value]
        except:
            new_value = 'N/A'
        new_dict[key] = new_value
    return new_dict

def valid_file( filepath ):
    if os.path.isfile( filepath ) and (os.path.getsize( filepath ) > 0):
        return True
    return False

def check_output_file( filepath ):
    if valid_file( filepath ):
        return
    else:
        msg = 'Expected output file not found! "{0}"'.format(filepath)
        log.error( msg )
        raise IOError( msg )

def remove_file( filepath ):
    if os.path.isfile( filepath ):
        try:
            os.remove( filepath )
        except:
            basename = os.path.basename( filepath )
            msg = 'Could not delete file! "%s"' % basename
            log.error( msg )
            raise IOError( msg )

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
