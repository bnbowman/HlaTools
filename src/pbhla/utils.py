import os
import subprocess
import logging
import string
import random

from collections import namedtuple

log = logging.getLogger()

BlasrM1 = namedtuple('BlasrM1', 'qname tname qstrand tstrand score pctsimilarity tstart tend tlength qstart qend qlength ncells')
BlasrM4 = namedtuple('BlasrM4', 'qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart tend tseqlength mapqv ncells clusterScore probscore numSigClusters')
BlasrM5 = namedtuple('BlasrM5', 'qname qlength z1 qalength qstrand tname tlength z2 talength tstrand score nmis nins ndel zscore qseq matchvector tseq')

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

def make_rand_string(minlength=6,maxlength=8):  
    length=random.randint(minlength,maxlength)  
    letters=string.ascii_letters+string.digits # alphanumeric, upper and lowercase  
    return ''.join([random.choice(letters) for _ in range(length)]) 

def getbash(cmd):
    out = subprocess.check_output(cmd, shell=True)
    return out  #This is the stdout from the shell command

def runbash(cmd):
    p = subprocess.call(cmd, 
                        shell=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
    return p

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

def check_output_file(filepath):
    try:
        assert os.path.isfile( filepath )
    except:
        msg = 'Expected output file not found! "{0}"'.format(filepath)
        log.info( msg )
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
