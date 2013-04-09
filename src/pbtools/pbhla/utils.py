import os
import string
import random
import subprocess
 
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

def create_directory( directory ):
    try:
        os.mkdir( directory )
    except:
        pass
