#! /usr/bin/env python
#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$

import os, re, sys
import logging
import subprocess

from pbcore.io.FastaIO import FastaReader

#from utils import which, create_directory

PULSE_METRICS = 'DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag'
NOISE_DATA = '-77.27,0.08654,0.00121'
PB_REGEX = 'm\d{6}_\d{6}_[a-zA-Z0-9]{4,6}_c\d{33}_s\d_p\d'

log = logging.getLogger()
log.setLevel( logging.INFO )

class ClusterResequencer(object):
    """
    A tool for resequencing clusters of rDNA sequences    
    """

    def __init__(self, read_file, 
                       ref_file, 
                       fofn_file, 
                       setup=None,
                       names=None,
                       output='Resequencing', 
                       nproc=1):
        self.read_file = read_file
        self.reference_file = ref_file
        self.fofn_file = fofn_file
        self._setup = setup
        self._names = names
        self._output = output
        self._nproc = nproc
        self._use_setup = None
        self.validate_settings()

    def validate_settings(self):
        # If the names are provided as a filename, parse it
        if isinstance(self._names, str):
            self._names = read_dict_file( self._names )
        if self._setup is not None:
            self._setup = os.path.abspath( self._setup )
        # Check that the fasta names can be mapped to PBI wells
        if (invalid_fasta_names( self.read_file ) and 
            invalid_dict_names( self._names )):
            msg = 'Cluster resequencing requires EITHER valid PacBio ' + \
                  'read-names OR a dictionary to convert them'
            raise Exception( msg )
        # Check for the availability of the various SMRT Analysis tools
        self.filter_plsh5 = which('filterPlsH5.py')
        self.compare_sequences = which('compareSequences.py')
        self.cmph5_tools = which('cmph5tools.py')
        self.load_pulses = which('loadPulses')
        self.variant_caller = which('variantCaller.py')
        # Check that either the tools or SMRT Analysis setup is available
        if all([ self.filter_plsh5,
                 self.compare_sequences,
                 self.cmph5_tools,
                 self.load_pulses,
                 self.variant_caller ]):
            self._use_setup = False
        elif self._setup and os.path.isfile( self._setup ):
            self._use_setup = True
            self._setup = os.path.abspath( self._setup )
            self.filter_plsh5 = 'filterPlsH5.py'
            self.compare_sequences = 'compareSequences.py'
            self.cmph5_tools = 'cmph5tools.py'
            self.load_pulses = 'loadPulses'
            self.variant_caller = 'variantCaller.py'
        else:
            msg = 'Cluster resequencing requires EITHER valid copies of ' + \
                  'the SMRT Analysis tools in the local path OR the ' + \
                  'path to a local SMRT Analysis setup script'
            raise Exception( msg )
        create_directory( self._output )

    def __call__(self):
        # Second we create a Rng.H5 file to mask other reads from Blasr
        whitelist_file = self.create_whitelist()
        rgnh5_fofn = self.create_rgnh5( whitelist_file )
        cmph5_file = self.create_cmph5( rgnh5_fofn )
        sorted_cmph5 = self.sort_cmph5( cmph5_file )
        self.run_load_pulses( sorted_cmph5 )
        #consensusFile = self.runQuiver( referenceFile, 
        #                                sortedCmpH5File,
        #                                count )

    def run_process(self, process_args, name):
        log.info("Executing child '%s' process" % name)
        if self._use_setup:
            process_args = ['source', self._setup, ';'] + process_args
        print process_args
        subprocess.call( process_args )

    def create_whitelist( self ):
        log.info('Creating cluster-specific whitelist file')
        whitelist_file = os.path.join(self._output, 'whitelist.txt')
        if os.path.exists( whitelist_file ):
            log.info('Existing whitelist file found, skipping...')
            return whitelist_file
        with open( whitelist_file, 'w' ) as handle:
            for record in FastaReader( self.read_file ):
                name = record.name.split()[0]
                if self._names:
                    name = self._names[name]
                handle.write(name + '\n')
        log.info('Finished writing the whitelist file')
        return os.path.abspath( whitelist_file )

    def create_rgnh5(self, whitelist_file):
        log.info('Creating cluster-specific RngH5 from Whitelist')
        output_dir = os.path.join(self._output, 'region_tables')
        output_fofn = os.path.join(self._output, 'region_tables.fofn')
        if os.path.exists( output_dir ) and os.path.exists( output_fofn ):
            log.info('Existing RngH5 detected, skipping...')
            return output_fofn
        process_args = [self.filter_plsh5, 
                        self.fofn_file,
                        '--outputDir=%s' % output_dir,
                        '--outputFofn=%s' % output_fofn,
                        '--filter=ReadWhitelist=%s' % whitelist_file]
        self.run_process( process_args, 'FilterPlsH5')
        log.info('Finished writing the cluster-specific RngH5')
        return output_fofn

    def update_fofn_paths(self, fofn_file):
        pass

    def create_cmph5(self, rgnh5_fofn ):
        log.info('Creating cluster-specific CmpH5')
        cmph5_file = os.path.join( self._output, 'cluster.cmp.h5' )
        if os.path.exists( cmph5_file ):
            log.info('Existing CmpH5 detected, skipping...')
            return cmph5_file
        process_args = [self.compare_sequences, 
                        '--useQuality',
                        '--h5pbi',
                        '--info',
                        '-x', '-bestn', '1',
                        '--nproc=%s' % self._nproc,
                        '--regionTable=%s' % rgnh5_fofn,
                        '--algorithm=blasr',
                        '--noiseData=%s' % NOISE_DATA,
                        '--h5fn=%s' % cmph5_file,
                        self.fofn_file,
                        self.reference_file]
        self.run_process( process_args, 'CompareSequences')
        log.info('Finished writing the cluster-specific CmpH5')
        return cmph5_file

    def sort_cmph5(self, cmph5_file ):
        log.info('Finished writing the cluster-specific CmpH5')
        sorted_cmph5 = os.path.join( self._output, 'cluster.sorted.cmp.h5')
        if os.path.exists( sorted_cmph5 ):
            log.info('Existing sorted CmpH5 detected, skipping...')
            return sorted_cmph5
        process_args = [self.cmph5_tools,
                        'sort',
                        '--outFile=%s' % sorted_cmph5,
                        cmph5_file]
        self.run_process( process_args, 'CmpH5Tools')
        log.info('Finished sorting the cluster-specific CmpH5')
        return sorted_cmph5

    def run_load_pulses(self, sorted_cmph5):
        log.info('Loading rich QV data into the CmpH5')
        process_args = [self.load_pulses,
                        self.fofn_file,
                        sorted_cmph5,
                        '-metrics', PULSE_METRICS]
        self.run_process( process_args, 'LoadPulses')
        log.info('Finished loading QV data into the CmpH5')

    def runQuiver(self, referenceFile, sortedCmpH5File, count):
        print "Running Quiver-consensus on Cluster #%s" % count
        consensusFile = 'cluster%s_consensus.fastq' % count
        if os.path.exists( consensusFile ):
            return consensusFile
        p = subprocess.Popen( [self.variantCaller,
                               '--algorithm=quiver',
                               '--numWorkers=%s' % self.numProc,
                               '--reference=%s' % referenceFile,
                               '--outputFile=%s' % consensusFile,
                               sortedCmpH5File])
        p.wait()
        return consensusFile

    def combineOutputSequences(self, outputSequenceFiles):
        print "Combining Consensus and Representative sequences"
        outputSequences = []
        for filename in outputSequenceFiles:
            for record in FastqReader( filename ):
                outputSequences.append( record )
        return outputSequences

    def outputCombinedSequences(self, combinedSequences ):
        print "Writing finished sequences to file"
        with FastqWriter( self.outputFile ) as handle:
            for record in combinedSequences:
                handle.writeRecord( record )

def invalid_fasta_names( fasta_file ):
    log.info('Checking read file for valid well names')
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

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    add = parser.add_argument
    add('read_file', metavar='READS', 
        help="Fasta or Fastq file of CCS sequences")
    add('ref_file',  metavar='REF',  
        help="Mothur list file of cluster data")
    add('fofn_file', metavar='FOFN', 
        help="BasH5 or FOFN of sequence data")
    add('--setup', metavar='SETUP_FILE',
        help='Path to the SMRT Analysis setup script')
    add('--names', metavar='DICT_FILE',
        help='Path to a dictionary file of sequence names')
    add('--output', default='reseq', metavar='DIR',
        help="Specify a directory for intermediate files")
    add('--nproc', type=int, default=1, metavar='INT',
        help="Number of processors to use")
    args = parser.parse_args()

    resequencer = ClusterResequencer( args.read_file, 
                                      args.ref_file, 
                                      args.fofn_file, 
                                      setup=args.setup,
                                      names=args.names,
                                      output=args.output, 
                                      nproc=args.nproc )
    resequencer()
