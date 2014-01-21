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
from pbhla.utils import which, create_directory
from pbhla.fasta.utils import invalid_fasta_names


PULSE_METRICS = 'DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag'
NOISE_DATA = '-77.27,0.08654,0.00121'
COVERAGE = 1000
CHEMISTRY = 'P4-C2.AllQVsMergingByChannelModel'

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
        self.read_file = read_file
        self.reference_file = ref_file
        self.fofn_file = fofn_file
        self._setup = setup
        self._output = output
        self._logs = os.path.join( self._output, 'logs' )
        self._scripts = os.path.join( self._output, 'scripts' )
        self._nproc = nproc
        self._use_setup = None
        self.validate_settings()

    def validate_settings(self):
        # If the names are provided as a filename, parse it
        if self._setup is not None:
            self._setup = os.path.abspath( self._setup )
        # Check that the fasta names can be mapped to PBI wells
        if invalid_fasta_names( self.read_file ):
            msg = 'Resequencer requires valid PacBio read-names'
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
            log.info('All required SMRT Analysis tools detected')
            self._use_setup = False
        elif self._setup and os.path.isfile( self._setup ):
            log.info('SMRT Analysis tools not detected, using Setup script')
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
        create_directory( self._scripts )
        create_directory( self._logs )

    def __call__(self):
        # Second we create a Rng.H5 file to mask other reads from Blasr
        whitelist_file = self.create_whitelist()
        rgnh5_fofn = self.create_rgnh5( whitelist_file )
        cmph5_file = self.create_cmph5( rgnh5_fofn )
        sorted_cmph5 = self.sort_cmph5( cmph5_file )
        self.run_load_pulses( cmph5_file, sorted_cmph5 )
        consensus_file = self.run_quiver( sorted_cmph5 )
        return consensus_file

    def run_process(self, process_args, name):
        log.info("Executing child '%s' process" % name)
        if self._use_setup:
            log.info('Executing subprocess indirectly via Shell Script')
            script = self.write_script( process_args, name)
            log_path = os.path.join( self._logs, name + '.log' )
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
        script_path = os.path.join( self._scripts, name + '_script.sh' )
        with open( script_path, 'w') as handle:
            handle.write('source %s\n' % self._setup)
            handle.write( ' '.join(process_args) + '\n' )
        return script_path

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
                        '--filter=ReadWhitelist=%s' % whitelist_file]
        self.run_process( process_args, 'FilterPlsH5')
        log.info('Finished writing the cluster-specific RngH5')
        return output_fofn

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

    def run_load_pulses(self, cmph5_file, sorted_cmph5):
        log.info('Loading rich QV data into the CmpH5')
        if os.path.getsize( sorted_cmph5 ) > 5*os.path.getsize( cmph5_file ):
            log.info('Existing CmpH5 appears to have QV data, skipping...')
            return
        process_args = [self.load_pulses,
                        self.fofn_file,
                        sorted_cmph5,
                        '-metrics', PULSE_METRICS]
        self.run_process( process_args, 'LoadPulses')
        log.info('Finished loading QV data into the CmpH5')

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
        help="Fasta or Fastq file of CCS sequences")
    add('ref_file',  metavar='REF',  
        help="Mothur list file of cluster data")
    add('fofn_file', metavar='FOFN', 
        help="BasH5 or FOFN of sequence data")
    add('--setup', metavar='SETUP_FILE',
        help='Path to the SMRT Analysis setup script')
    add('--output', default='reseq', metavar='DIR',
        help="Specify a directory for intermediate files")
    add('--nproc', type=int, default=1, metavar='INT',
        help="Number of processors to use")

    args = parser.parse_args()
    
    # If the args were valid, initialize a logger before running
    logging.basicConfig( level=logging.INFO )

    # Run the specified resequencing process
    resequencer = ClusterResequencer( args.read_file, 
                                      args.ref_file, 
                                      args.fofn_file, 
                                      setup=args.setup,
                                      output=args.output, 
                                      nproc=args.nproc )
    resequencer()
