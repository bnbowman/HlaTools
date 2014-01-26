#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import logging, logging.config
import subprocess

from pbhla import __LOG__
from pbhla.utils import which

logging.config.fileConfig( __LOG__ )
log = logging.getLogger(__name__)

class ExternalTool( object ):
    """An abstract class for running external command-line tools"""

    def __init__(self, setup=None):
        self._setup = setup
        self._use_setup = self.test_setup()
        self._args

    def test_setup(self):
        """Determine whether we need a setup script, and which environment to use"""
        if self._setup and which(self.exe):
            # TODO: Add validation that the setup script works
            log.info('"%s" not detected, using supplied environment' % self.name)
            return True
        elif self.setup and which(self.exe) is None:
            log.info('%s detected, but using supplied environment instead.' % self.name + \
                     'Do not pass a "setup" argument to use local environment')
            return True
        elif self.setup is None and which(self.exe) is None:
            log.info('"%s" detected, using local environment' % self.name)
            return False
        else:
            msg = '"%s" requires EITHER a valid executable in the local ' % self.name + \
                  'path OR a virtualenv setup script'
            log.error( msg )
            raise Exception( msg )

    @property
    def name(self):
        raise NotImplementedError("Subclasses should implement this!")

    @property
    def exe(self):
        raise NotImplementedError("Subclasses should implement this!")

    @property
    def setup(self):
        return self._setup

    @property
    def use_setup(self):
        return self._use_setup

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
