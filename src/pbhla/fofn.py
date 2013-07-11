import os, logging

from pbhla.utils import write_list_file

log = logging.getLogger()

def create_baxh5_fofn( input_file, output_file ):
    if input_file.endswith('.fofn'):
        baxh5_files = _parse_fofn( input_file )
    elif input_file.endswith('.bas.h5'):
        baxh5_files = _parse_bash5( input_file )
    elif input_file.endswith('.bax.h5'):
        baxh5_files = [input_file]
    elif input_file.endswith('.fa') or input_file.endswith('.fasta'):
        baxh5_files = []
    else:
        msg = 'Invalid input filetype "%s"' % input_file
        log.info( msg )
        raise TypeError( msg )
    write_list_file( baxh5_files, output_file )
    return output_file

def _parse_fofn( filename ):
    baxh5_files = []
    with open( filename, 'r' ) as handle:
        for line in handle:
            filename = line.strip()
            if not filename:
                continue
            if filename.endswith('.bas.h5'):
                baxh5_files += _parse_bash5( filename )
            if filename.endswith('.bax.h5'):
                baxh5_files.append( filename )
    return baxh5_files

def _parse_bash5( filename ):
    root = '.'.join(filename.split('.')[:-2])
    bax_names = ['%s.%s.bax.h5' % (root, i) for i in range(1,4)]
    all_bax = all([os.path.isfile(bax) for bax in bax_names])
    any_bax = any([os.path.isfile(bax) for bax in bax_names])
    if all_bax:
        return bax_names
    elif not any_bax:
        return [filename]
    else:
        msg = 'Mixture of valid and invalid Bax.h5 files!'
        log.info( msg )
        raise ValueError( msg )
