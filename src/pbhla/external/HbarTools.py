#! /usr/bin/env python
import logging, subprocess, os

from shutil import copy

MIN_LENGTH = 2500
MIN_SCORE = 0.80

class HbarRunner( object ):

    def __init__(self, input_file, output_dir, config_file=None, 
                                               min_length=None, 
                                               min_score=None):
        self.input_file = input_file
        self.origin_dir = os.getcwd()
        self.output_dir = os.path.abspath( output_dir )
        self.config_file = config_file
        self.log_file = input_file.split('.')[0] + '.log'
        self.min_length = min_length or MIN_LENGTH
        self.min_score = min_score or MIN_SCORE
        self.initialize_logging()

    def initialize_logging(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing HbarRunner with the following:")
        self.log.info("\tInput: {0}".format(self.input_file))
        self.log.info("\tOutput: {0}".format(self.output_dir))
        self.log.info("\tConfig: {0}".format(self.config_file))
        self.log.info("\tLog: {0}".format(self.log_file))
        self.log.info("\tMinimum Length: {0}".format(self.min_length))
        self.log.info("\tMinimum Score: {0}".format(self.min_score))

    def __call__(self):
        os.chdir( self.output_dir )
        self.write_config()
        self.run_hbar()
        os.chdir( self.origin_dir )

    def write_config(self):
        self.config_output = os.path.join( self.output_dir, 'HLA_HBAR.cfg' )
        self.log.info("Writing HBAR configuration to {0}".format(self.config_output))
        if self.config_file:
            self.log.info("Copying existing configuration from {0}".format(self.config_file))
            shutil.copy( self.config_file, self.config_output )
        else:
            self.log.info("Creating new configuration from options")
            with open(self.config_output, 'w') as handle:
                handle.write( self.hgap_config )

    def run_hbar(self):
        self.log.info("Running HBAR Workflow")
        hbar_args = ["HBAR_WF2.py", self.config_output]
        with open(self.log_file, 'w') as log_handle:
            subprocess.check_call( hbar_args, stdout=log_handle )
        self.log.info("HBAR completed successfully")

    @property
    def hgap_config(self):
        return "\n".join(['[General]',
            '#list of files of the initial bas.h5 files',
            'input_fofn = {0}'.format(self.input_file),
            '#The length cutoff used for seed reads',
            'length_cutoff = {0}'.format(MIN_LENGTH),
            '#The length cutoff used for pre-assembly',
            'length_cutoff_pr = {0}'.format(MIN_LENGTH),
            '#The read quality cutoff used for seed reads',
            'RQ_threshold = {0}'.format(MIN_SCORE),
            '#SGE job option for distributed mapping',
            'sge_option_dm = -pe smp 12 -q secondary',
            '#SGE job option for pre-assembly',
            'sge_option_qf = -pe smp 1 -q secondary',
            '#SGE job option for pre-assembly',
            'sge_option_pa = -pe smp 16 -q secondary',
            '#SGE job option for CA',
            'sge_option_ca = -pe smp 4 -q secondary',
            '#SGE job option for Quiver',
            'sge_option_qv = -pe smp 16 -q secondary',
            '#SGE job option for "qsub -sync y" to sync stage',
            'sge_option_ck = -pe smp 1 -q secondary',
            '# Options for blasr mapping',
            'blasr_opt = -nCandidates 32 -minMatch 12 -maxLCPLength 15 -bestn 16 -minPctIdentity 75.0 -maxScore -1000 -nproc 16',
            '#This is used for running quiver',
            'SEYMOUR_HOME = /mnt/secondary/Smrtpipe/builds/Assembly_Mainline_Nightly_Archive/build470-116466/',
            '#The number of best alignment hits used for pre-assembly',
            '#It should be about the same as the final PLR coverage, slight higher might be OK.',
            'bestn = 32',
            '# target choices are "pre_assembly", "draft_assembly", "all"',
            '# "pre_assembly" : generate pre_assembly for any long read assembler to use',
            '# "draft_assembly": automatic submit CA assembly job when pre-assembly is done',
            '# "all" : submit job for using Quiver to do final polish',
            'target = draft_assembly',
            '# number of chunks for distributed mapping',
            'preassembly_num_chunk = 8',
            '# number of chunks for pre-assembly',
            'dist_map_num_chunk = 6',
            '# Multi-file chunking parameters',
            'q_chunk_size = 1',
            't_chunk_size = 1',
            '# "tmpdir" is for preassembly, better in a ramdisk',
            'tmpdir = /tmp',
            '# "big_tmpdir" is for quiver, better in a big disk',
            'big_tmpdir = /tmp',
            '# Other parameters',
            'min_cov = 8',
            'max_cov = 64',
            'trim_align = 1',
            'trim_plr = 1',
            'q_nproc = 16'])
