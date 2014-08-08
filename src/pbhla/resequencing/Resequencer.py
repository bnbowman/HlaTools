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
import logging, logging.config
import subprocess

from pbcore.io.FastaIO import FastaReader
from pbhla.fasta.utils import is_pacbio_record
from pbhla.utils import which, create_directory, validate_file

log = logging.getLogger()

PULSE_METRICS = 'DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag'
COVERAGE = 200
CHEMISTRY = 'P4-C2.AllQVsMergingByChannelModel'

class Resequencer(object):
    """
    A tool for automatic the resequencing small subsets of a PacBio dataset
    """

    def __init__(self, setup=None, nproc=1, debug=False):
        """Initialize cross-cluster and object-specific settings"""
        self._setup = self.validate_setup_path( setup ) if setup is not None else None
        self._nproc = nproc
        self._debug = debug
        if debug:
            log.setLevel( logging.DEBUG )
        log.debug("TESTING")

        # Find the various required scripts and determine if
        self.samtools       = self.validate_program( ['samtools'] )
        self.filter_plsh5   = self.validate_program( ['filter_plsh5.py', 'filterPlsH5.py'] )
        self.pbalign        = self.validate_program( ['pbalign'] )
        self.variant_caller = self.validate_program( ['variantCaller.py'] )
        self._use_setup     = self.validate_setup()

        log.info( self.filter_plsh5 )
        log.info( self.pbalign )
        log.info( self.variant_caller )

        # Initialize properties to be used for each call
        self._counter = 0
        self._output = None
        self._logs = None
        self._scripts = None

    @property
    def setup(self):
        return self._setup

    @property
    def nproc(self):
        return self._nproc

    @property
    def debug(self):
        return self._debug

    @property
    def subprograms_in_path(self):
        if self.samtools.startswith( self.setup ):
            return False
        if self.filter_plsh5.startswith( self.setup ):
            return False
        if self.pbalign.startswith( self.setup ):
            return False
        if self.variant_caller.startswith( self.setup ):
            return False
        return True

    def validate_setup_path(self, setup):
        setup = os.path.abspath( setup )
        if setup.endswith('/etc/setup.sh'):
            setup = setup[:-13]
            print setup
        setup_script = os.path.join( setup, 'etc/setup.sh')
        print setup_script
        if not validate_file( setup_script ):
            msg = "Supplied setup-script does not appear to be valid"
            log.error( msg )
            raise ValueError( msg )
        return setup

    def validate_program( self, programs):
        # First we try to find the program in the local path
        program_paths = [which(program) for program in programs]
        program_paths = [path for path in program_paths if path is not None]
        print programs
        print program_paths
        if program_paths:
            return program_paths[0]
        elif self.setup is None:
            msg = 'No program from the set "%s" found in PATH and no SMRT Analysis supplied' % str(programs)
            log.error( msg )
            raise IOError( msg )

        # Fallback to the Setup Script if needed
        setup_programs = [os.path.join( self.setup, 'analysis/bin', p ) for p in programs]
        setup_programs = [validate_file( p ) for p in setup_programs]
        setup_programs = [path for path in setup_programs if path is not None]
        if setup_programs:
            return setup_programs[0]
        else:
            msg = 'No program from the set "%s" found in either PATH or SMRT Analysis env' % str(programs)
            log.error( msg )
            raise IOError( msg )

    def validate_setup(self):
        """Determine whether we need a setup script, and which environment to use"""
        if self.subprograms_in_path:
            need_setup = False
        else:
            need_setup = True

        # Determine whether we can proceed or not
        if self._setup and need_setup:
            # TODO: Add validation that the setup script works
            log.info('SMRT Analysis tools not detected, using supplied '
                     'SMRT Analysis environment')
            return True
        elif self._setup and not need_setup:
            log.info('SMRT Analysis tools were detected, but using supplied '
                     'SMRT Analysis environment instead.  Do not pass a '
                     '"setup" argument to use local environment')
            return True
        elif self._setup is None and not need_setup:
            log.info('SMRT Analysis tools detected, using local '
                     'SMRT Analysis environment')
            return False
        else:
            msg = 'Cluster resequencing requires EITHER valid copies of ' + \
                  'SMRT Analysis tools in the local path OR the ' + \
                  'path to a local SMRT Analysis setup script'
            log.error( msg )
            raise Exception( msg )


    def initialize_output(self, output ):
        """Initialize the cluster-specific output folders"""
        # TODO: Check for existing directories and do something
        self._output = os.path.abspath( output )
        create_directory( self._output )
        self._scripts = os.path.join( self._output, 'scripts' )
        create_directory( self._scripts )
        self._logs = os.path.join( self._output, 'logs' )
        create_directory( self._logs )


    def run_process(self, process_args, name):
        log.info("Executing child '%s' process" % name)
        if self._use_setup:
            log.info('Executing subprocess indirectly via Shell Script')
            script = self.write_script( process_args, name )
            log_path = self.get_log_path( name )
            with open( log_path, 'w' ) as log_handle:
                p = subprocess.Popen( ['source', script], 
                                       executable='/bin/bash',
                                       stderr=subprocess.STDOUT,
                                       stdout=log_handle)
                p.wait()
        else:
            log.info('Executing subprocess directly via Subprocess')
            p = subprocess.Popen( process_args )
            p.wait()
        log.info('Child process finished successfully')


    def write_script( self, process_args, name ):
        script_path = self.get_script_path( name )
        with open( script_path, 'w') as handle:
            handle.write('export SEYMOUR_HOME=%s\n' % self._setup)
            handle.write('source %s/etc/setup.sh\n' % self._setup)
            arg_string = ' '.join([str(arg) for arg in process_args])
            handle.write( arg_string + '\n' )
        return script_path


    def get_script_path( self, name ):
        self._counter += 1
        script_name = '%s_%s_script.sh' % (self._counter, name )
        return os.path.join( self._scripts, script_name )


    def get_log_path( self, name ):
        log_name = '%s_%s.log' % (self._counter, name)
        return os.path.join( self._logs, log_name )


    def create_whitelist( self, read_file ):
        log.info('Creating cluster-specific whitelist file')
        if read_file.endswith('.txt') or read_file.endswith('.zmws'):
            log.info('Whitelist is already in a valid format, skipping...')
            return read_file
        whitelist_file = os.path.join(self._output, 'whitelist.txt')
        if os.path.exists( whitelist_file ):
            log.info('Existing whitelist file found, skipping...')
            return whitelist_file

        # Parse the individual ZMWs of interest
        zmws = set()
        for record in FastaReader( read_file ):
            assert is_pacbio_record( record )
            name = record.name.split()[0]
            zmw = '/'.join( name.split('/')[:2] )
            zmws.add( zmw )

        # Write the ZMWs to file
        with open( whitelist_file, 'w' ) as handle:
            for zmw in sorted( zmws ):
                handle.write(zmw + '\n')
        log.info('Finished writing the whitelist file\n')
        return os.path.abspath( whitelist_file )

    def create_fai( self, reference_file ):
        log.info('Creating reference-specific index file')
        index_file = reference_file + '.fai'
        if os.path.exists( index_file ):
            log.info('Existing reference detected, skipping...')
            self._counter += 1
            return index_file
        process_args = [self.samtools,
                        'faidx',
                        reference_file]
        self.run_process( process_args, 'SamTools' )
        log.info('Finished writing the the reference index file')
        return index_file

    def create_rgnh5( self, data_file, whitelist_file ):
        log.info('Creating cluster-specific RngH5 from Whitelist')
        output_dir = os.path.join(self._output, 'region_tables')
        output_fofn = os.path.join(self._output, 'region_tables.fofn')
        if os.path.exists( output_dir ) and os.path.exists( output_fofn ):
            log.info('Existing RngH5 detected, skipping...')
            self._counter += 1
            return output_fofn
        process_args = [self.filter_plsh5, 
                        data_file,
                        '--outputDir=%s' % output_dir,
                        '--outputFofn=%s' % output_fofn,
                        '--filter="ReadWhitelist=%s,MinReadScore=0.75,MaxSRL=3600"' % whitelist_file]
        self.run_process( process_args, 'FilterPlsH5' )
        log.info('Finished writing the cluster-specific RngH5')
        return output_fofn

    def run_pbalign( self, data_file, reference_file, rgnh5_file, min_length ):
        log.info('Creating cluster-specific CmpH5')
        cmph5_file = os.path.join( self._output, 'cluster.cmp.h5' )
        if os.path.exists( cmph5_file ):
            log.info('Existing CmpH5 detected, skipping...')
            self._counter += 1
            return cmph5_file
        process_args = [self.pbalign,
                        data_file,
                        reference_file,
                        cmph5_file,
                        '--regionTable=%s' % rgnh5_file,
                        '--forQuiver',
                        '--minLength=%s' % min_length,
                        '--hitPolicy=allbest',
                        '--nproc=%s' % self._nproc]
        self.run_process( process_args, 'PBAlign')
        log.info('Finished writing the cluster-specific CmpH5')
        return cmph5_file


    def run_quiver( self, cmph5_file, reference_file ):
        log.info('Running Quiver to generate HQ consensus')
        consensus_fastq = os.path.join( self._output, 'consensus.fastq')
        consensus_fasta = os.path.join( self._output, 'consensus.fasta')
        if ( os.path.exists( consensus_fastq ) and
             os.path.isfile( consensus_fasta )):
            log.info('Existing consensus found, skipping...')
            self._counter += 1
            return consensus_fastq, consensus_fasta
        process_args = [self.variant_caller,
                        '--algorithm=quiver',
                        '--verbose',
                        #'--parameterSet=%s' % CHEMISTRY,
                        '--numWorkers=%s' % self._nproc,
                        '--reference=%s' % reference_file,
                        '--coverage=%s' % COVERAGE,
                        '--outputFile=%s' % consensus_fastq,
                        '--outputFile=%s' % consensus_fasta,
                        cmph5_file]
        self.run_process( process_args, 'Quiver' )
        log.info('Finished creating the Quiver consensus')
        return consensus_fastq, consensus_fasta


    def __call__(self, data_file, read_file, reference_file, output='Resequencing', min_length=500):
        # Validate and set-up run-specific values
        self._counter  = 0
        data_file      = validate_file( data_file )
        read_file      = validate_file( read_file )
        reference_file = validate_file( reference_file )
        self.initialize_output( output )

        # Next we run the 5-step Quiver Process
        whitelist_file = self.create_whitelist( read_file )
        self.create_fai( reference_file )
        rgnh5_file = self.create_rgnh5( data_file, whitelist_file )
        cmph5_file = self.run_pbalign( data_file, reference_file, rgnh5_file, min_length )
        consensus_file = self.run_quiver( cmph5_file, reference_file )
        return consensus_file


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    add = parser.add_argument
    add('read_file', metavar='READS', 
        help="Fasta file of Cluster-specific reads")
    add('ref_file',  metavar='REF',  
        help="Fasta file of sequences to polish with Quiver")
    add('fofn_file', metavar='FOFN', 
        help="BasH5 or FOFN of sequence data")
    add('--setup', metavar='SETUP_FILE',
        help='Path to the SMRT Analysis setup script')
    add('--output', default='resequencing', metavar='DIR',
        help="Specify a directory for intermediate files")
    add('--nproc', type=int, default=1, metavar='INT',
        help="Number of processors to use")
    add('--debug', action='store_true',
        help="Flag to enable Debug mode")
    args = parser.parse_args()

    # Run the specified resequencing process
    resequencer = Resequencer( setup=args.setup,
                               nproc=args.nproc,
                               debug=args.debug )
    resequencer( args.fofn_file,
                 args.read_file,
                 args.ref_file,
                 output=args.output )
