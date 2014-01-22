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
from pbhla import __LOG__
from pbhla.fasta.utils import invalid_fasta_names
from pbhla.utils import which, create_directory

PULSE_METRICS = 'DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag'
NOISE_DATA = '-77.27,0.08654,0.00121'
COVERAGE = 200
CHEMISTRY = 'P4-C2.AllQVsMergingByChannelModel'

logging.config.fileConfig( __LOG__ )
log = logging.getLogger()

class ClusterResequencer(object):
    """
    A tool for resequencing clusters of rDNA sequences    
    """

    def __init__(self, read_file, 
                       ref_file, 
                       fofn_file, 
                       setup=None,
                       output='Resequencing', 
                       nproc=1):
        self._read_file = read_file
        self._reference_file = ref_file
        self._fofn_file = fofn_file
        self._counter = 0
        self._setup = os.path.abspath( setup ) if setup is not None else None
        self._output = output
        self._logs = os.path.join( self._output, 'logs' )
        self._scripts = os.path.join( self._output, 'scripts' )
        self._nproc = nproc
        self._use_setup = None
        self.validate_settings()

    def validate_settings(self):
        # Check that the fasta names can be mapped to PBI wells
        if invalid_fasta_names( self.read_file ):
            msg = 'Resequencer requires valid PacBio read-names'
            log.error( msg )
            raise ValueError( msg )

        # Check for the availability of the various SMRT Analysis tools
        self.filter_plsh5 = which('filterPlsH5.py')
        self.pbalign = which('pbalign.py')
        self.variant_caller = which('variantCaller.py')

        # Check that either the tools or SMRT Analysis setup is available
        if all([ self.filter_plsh5,
                 self.pbalign,
                 self.variant_caller ]):
            log.info('All required SMRT Analysis tools detected')
            self._use_setup = False
        elif self._setup and os.path.isfile( self._setup ):
            log.info('SMRT Analysis tools not detected, using Setup script')
            self._use_setup = True
            self._setup = os.path.abspath( self._setup )
            self.filter_plsh5 = 'filterPlsH5.py'
            self.pbalign = 'pbalign.py'
            self.variant_caller = 'variantCaller.py'
        else:
            msg = 'Cluster resequencing requires EITHER valid copies of ' + \
                  'the SMRT Analysis tools in the local path OR the ' + \
                  'path to a local SMRT Analysis setup script'
            log.error( msg )
            raise Exception( msg )

        create_directory( self._output )
        create_directory( self._scripts )
        create_directory( self._logs )

    @property
    def read_file(self):
        return self._read_file

    @property
    def reference_file(self):
        return self._reference_file

    @property
    def fofn_file(self):
        return self._fofn_file

    def __call__(self):
        # Second we create a Rng.H5 file to mask other reads from Blasr
        whitelist_file = self.create_whitelist()
        rgnh5_fofn = self.create_rgnh5( whitelist_file )
        cmph5_file = self.run_pbalign( rgnh5_fofn )
        consensus_file = self.run_quiver( cmph5_file )
        return consensus_file

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
            handle.write('source %s\n' % self._setup)
            handle.write( ' '.join(process_args) + '\n' )
        return script_path

    def get_script_path( self, name ):
        self._counter += 1
        script_name = '%s_%s_script.sh' % (self._counter, name )
        return os.path.join( self._scripts, script_name )

    def get_log_path( self, name ):
        log_name = '%s_%s.log' % (self._counter, name)
        return os.path.join( self._logs, log_name )

    def create_whitelist( self ):
        log.info('Creating cluster-specific whitelist file')
        whitelist_file = os.path.join(self._output, 'whitelist.txt')
        if os.path.exists( whitelist_file ):
            log.info('Existing whitelist file found, skipping...')
            return whitelist_file
        # Parse the individual ZMWs of interest
        zmws = []
        for record in FastaReader( self.read_file ):
            name = record.name.split()[0]
            zmw = '/'.join( name.split('/')[:2] )
            zmws.append( zmw )
        # Write the ZMWs to file
        with open( whitelist_file, 'w' ) as handle:
            for zmw in sorted(set( zmws )):
                handle.write(zmw + '\n')
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
                        '--filter="ReadWhitelist=%s,MinSRL=1500,MinReadScore=0.8,MaxSRL=3700"' % whitelist_file]
        self.run_process( process_args, 'FilterPlsH5')
        log.info('Finished writing the cluster-specific RngH5')
        return output_fofn

    def run_pbalign(self, rgnh5_fofn ):
        log.info('Creating cluster-specific CmpH5')
        cmph5_file = os.path.join( self._output, 'cluster.cmp.h5' )
        if os.path.exists( cmph5_file ):
            log.info('Existing CmpH5 detected, skipping...')
            return cmph5_file
        process_args = [self.pbalign,
                        self.fofn_file,
                        self.reference_file,
                        cmph5_file,
                        '--regionTable=%s' % rgnh5_fofn,
                        '--forQuiver',
                        '--hitPolicy=randombest',
                        '--nproc=%s' % self._nproc]
        self.run_process( process_args, 'CompareSequences')
        log.info('Finished writing the cluster-specific CmpH5')
        return cmph5_file

    def run_quiver(self, sorted_cmph5):
        log.info('Running Quiver to generate HQ consensus')
        consensus_fastq = os.path.join( self._output, 'consensus.fastq')
        consensus_fasta = os.path.join( self._output, 'consensus.fasta')
        if ( os.path.exists( consensus_fastq ) and
             os.path.isfile( consensus_fasta )):
            log.info('Existing consensus found, skipping...')
            return consensus_fastq, consensus_fasta
        process_args = [self.variant_caller,
                        '--algorithm=quiver',
                        '--verbose',
                        '--parameterSet=%s' % CHEMISTRY,
                        '--numWorkers=%s' % self._nproc,
                        '--reference=%s' % self.reference_file,
                        '--coverage=%s' % COVERAGE,
                        '--outputFile=%s' % consensus_fastq,
                        '--outputFile=%s' % consensus_fasta,
                        sorted_cmph5]
        self.run_process( process_args, 'Quiver' )
        log.info('Finished creating the Quiver consensus')
        return consensus_fastq, consensus_fasta



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
    args = parser.parse_args()

    # Run the specified resequencing process
    resequencer = ClusterResequencer( args.read_file, 
                                      args.ref_file, 
                                      args.fofn_file, 
                                      setup=args.setup,
                                      output=args.output, 
                                      nproc=args.nproc )
    resequencer()
