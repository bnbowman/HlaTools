import os
import subprocess
import logging
import string
import random

from collections import namedtuple

log = logging.getLogger()

BlasrM1 = namedtuple('BlasrM1', 'qname tname qstrand tstrand score pctsimilarity tstart tend tlength qstart qend qlength ncells')

def write_fofn( file_list, output_file ):
    with open(output_file, 'w') as handle:
        for filename in file_list:
            print >> handle, filename

def read_fofn( fofn_file ):
    file_list = []
    with open(fofn_file, 'r') as handle:
        for line in handle:
            file_list.append( line.strip().split()[0] )
    return file_list

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
        try:
            new_value = ref_dict[old_value]
        except:
            new_value = 'N/A'
        new_dict[key] = new_value
    return new_dict

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
