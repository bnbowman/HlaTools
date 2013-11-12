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

import os
import logging
import subprocess

from pbcore.io.FastaIO import FastaReader
from pbhla.log import initialize_logger
from pbhla.utils import create_directory, valid_file, check_output_file

PULSE_METRICS = 'DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag'
NOISE_DATA = '-77.27,0.08654,0.00121'
CHEMISTRY = 'P4-C2.AllQVsMergingByChannelModel'

log = logging.getLogger()

class Resequencer(object):
    """
    A tool for resequencing clusters of rDNA sequences    
    """

    def __init__(self, fofn_file, setup=None, nproc=1):
        self.fofn_file = fofn_file
        self._setup = setup
        self._nproc = nproc
        self._use_setup = None
        self._tools = {}
        self._validate_settings()

    def __getattr__(self, item):
        return self._tools[item]

    def _validate_settings(self):
        """
        Confirm that all of the required utilities are accessible,
            or barring that, that a valid Setup script was supplied
        """
        # Make sure we have the absolute path to the setup
        if self._setup is not None:
            self._setup = os.path.abspath( self._setup )
        # Check for the availability of the various SMRT Analysis tools
        self._tools['filter_plsh5'] = which('filterPlsH5.py')
        self._tools['compare_sequences'] = which('compareSequences.py')
        self._tools['cmph5_tools'] = which('cmph5tools.py')
        self._tools['load_pulses'] = which('loadPulses')
        self._tools['variant_caller'] = which('variantCaller.py')
        # Check that either the tools or SMRT Analysis setup is available
        if all([ self.filter_plsh5,
                 self.compare_sequences,
                 self.cmph5_tools,
                 self.load_pulses,
                 self.variant_caller ]):
            log.info('All required SMRT Analysis tools detected')
            self._use_setup = False
        elif self._setup and os.path.isfile( self._setup ):
            log.info('Some SMRT Analysis tools not detected, using Setup script')
            self._use_setup = True
            self._tools['filter_plsh5'] = 'filterPlsH5.py'
            self._tools['compare_sequences'] = 'compareSequences.py'
            self._tools['cmph5_tools'] = 'cmph5tools.py'
            self._tools['load_pulses'] = 'loadPulses'
            self._tools['variant_caller'] = 'variantCaller.py'
        else:
            msg = 'Resequencing requires EITHER valid copies of ' + \
                  'the SMRT Analysis tools in the local path OR the ' + \
                  'path to a local SMRT Analysis setup script'
            log.error( msg )
            raise Exception( msg )

    def __call__( self, subread_file, reference_file, output=None ):
        """
        Execute all the necessary steps to Quiver-correct a reference file
            with the subreads from a different file
        """
        # Setup the output directory
        output = output or _get_output( subread_file )
        create_directory( output )
        # Run the 6 parts of the Quiver resequencing work-flow
        whitelist = self.create_whitelist( subread_file, output )
        rgnh5_fofn = self.create_rgnh5( whitelist, output )
        cmph5_file = self.create_cmph5( rgnh5_fofn, reference_file, output )
        sorted_cmph5 = self.sort_cmph5( cmph5_file, output )
        self.run_load_pulses( cmph5_file, sorted_cmph5, output )
        return self.run_quiver( sorted_cmph5, reference_file, output )

    def run_process( self, process_args, name, output ):
        """
        Execute a sub-process via Popen
        """
        log_path = os.path.join( output, name + '.log' )
        if self._use_setup:
            log.info('Executing %s subprocess indirectly via shell-script' % name)
            script = self.write_script( process_args, name, output )
            with open( log_path, 'w' ) as log_handle:
                p = subprocess.Popen( ['source', script], 
                                       executable='/bin/bash',
                                       stderr=subprocess.STDOUT,
                                       stdout=log_handle )
                p.wait()
        else:
            log.info('Executing %s subprocess directly' % name)
            with open( log_path, 'w' ) as log_handle:
                p = subprocess.Popen( process_args,
                                      stderr=subprocess.STDOUT,
                                      stdout=log_handle )
                p.wait()
        log.info('%s subprocess finished successfully' % name)

    def write_script( self, process_args, name, output ):
        """
        Convert a supplied argument list into a Bash script
        """
        script_path = os.path.join( output, name + '_script.sh' )
        with open( script_path, 'w') as handle:
            handle.write('source %s\n' % self._setup)
            handle.write( ' '.join(process_args) + '\n' )
        return script_path

    def create_whitelist( self, subread_file, output ):
        log.info('Creating target-specific whitelist file')
        whitelist_file = os.path.join( output, 'whitelist.txt')
        if valid_file( whitelist_file ):
            log.info('Existing whitelist file found, skipping...')
            return whitelist_file
        with open( whitelist_file, 'w' ) as handle:
            zmw_set = set()
            for record in FastaReader( subread_file ):
                name = record.name.split()[0]
                zmw = '/'.join( name.split('/')[:2] )
                if zmw in zmw_set:
                    continue
                else:
                    zmw_set.add( zmw )
                    handle.write( zmw + '\n' )
        check_output_file( whitelist_file )
        log.info('Finished writing the whitelist file')
        return whitelist_file

    def create_rgnh5(self, whitelist, output):
        """
        Create a RngH5 corresponding to the supplied subreads
        """
        log.info('Creating target-specific RngH5 from Whitelist')
        output_dir = os.path.join( output, 'region_tables' )
        output_fofn = os.path.join( output, 'region_tables.fofn' )
        if valid_file( output_fofn ):
            log.info('Existing RngH5 detected, skipping...')
            return output_fofn
        process_args = [self.filter_plsh5, 
                        self.fofn_file,
                        '--outputDir=%s' % output_dir,
                        '--outputFofn=%s' % output_fofn,
                        '--filter=ReadWhitelist=%s' % whitelist]
        self.run_process( process_args, 'FilterPlsH5', output )
        check_output_file( output_fofn )
        log.info('Finished writing the cluster-specific RngH5')
        return output_fofn

    def create_cmph5(self, rgnh5_fofn, reference_file, output ):
        """
        Create a CmpH5 using using only the supplied subreads/RngH5
        """
        log.info('Creating target-specific CmpH5')
        cmph5_file = os.path.join( output, 'locus.cmp.h5' )
        if valid_file( cmph5_file ):
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
                        reference_file]
        self.run_process( process_args, 'CompareSequences', output )
        check_output_file( cmph5_file )
        log.info('Finished writing the cluster-specific CmpH5')
        return cmph5_file

    def sort_cmph5(self, cmph5_file, output ):
        """
        Sort the CmpH5 file as a pre-requisite for Quiver
        """
        log.info('Creating a sorted target-specific CmpH5')
        sorted_cmph5 = os.path.join( output, 'locus.sorted.cmp.h5')
        if os.path.exists( sorted_cmph5 ):
            log.info('Existing sorted CmpH5 detected, skipping...')
            return sorted_cmph5
        process_args = [self.cmph5_tools,
                        'sort',
                        '--outFile=%s' % sorted_cmph5,
                        cmph5_file]
        self.run_process( process_args, 'CmpH5Tools', output )
        check_output_file( sorted_cmph5 )
        log.info('Finished sorting the cluster-specific CmpH5')
        return sorted_cmph5

    def run_load_pulses(self, cmph5_file, sorted_cmph5, output ):
        """
        Load the rich QV data into the CmpH5, also a quiver pre-req
        """
        log.info('Loading rich QV data into the CmpH5')
        if os.path.getsize( sorted_cmph5 ) > 5*os.path.getsize( cmph5_file ):
            log.info('Existing CmpH5 appears to have QV data, skipping...')
            return
        process_args = [self.load_pulses,
                        self.fofn_file,
                        sorted_cmph5,
                        '-metrics', PULSE_METRICS]
        self.run_process( process_args, 'LoadPulses', output )
        log.info('Finished loading QV data into the CmpH5')

    def run_quiver(self, sorted_cmph5, reference_file, output ):
        """
        Run Quiver on the prepared, target-specific CmpH5 file
        """
        log.info('Running Quiver to generate HQ consensus')
        consensus_fastq = os.path.join( output, 'consensus.fastq')
        consensus_fasta = os.path.join( output, 'consensus.fasta')
        if valid_file( consensus_fastq ):
            log.info('Existing consensus found, skipping...')
            return consensus_fastq
        process_args = [self.variant_caller,
                        '--algorithm=quiver',
                        '--verbose',
                        '--parameterSet=%s' % CHEMISTRY,
                        '--numWorkers=%s' % self._nproc,
                        '--reference=%s' % reference_file,
                        '--outputFile=%s' % consensus_fastq,
                        '--outputFile=%s' % consensus_fasta,
                        sorted_cmph5]
        self.run_process( process_args, 'Quiver', output )
        log.info('Finished creating the Quiver consensus')
        return consensus_fastq

def _get_output( subread_file ):
    return '.'.join( subread_file.split('.')[:-1] )

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

if __name__ == '__main__':
    import sys

    fofn_file = sys.argv[1]
    setup_file = sys.argv[2]
    subread_file = sys.argv[3]
    reference_file = sys.argv[4]

    # Run the specified resequencing process
    initialize_logger()
    resequencer = Resequencer( fofn_file, setup=setup_file, nproc=12 )
    resequencer( subread_file, reference_file )