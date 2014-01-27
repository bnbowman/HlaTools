#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import logging, logging.config
import subprocess
import tempfile
import shlex

from pbhla import __LOG__
from pbhla.utils import which

logging.config.fileConfig( __LOG__ )
log = logging.getLogger(__name__)

class ExternalTool( object ):
    """An abstract class for running external command-line tools"""
    counter = 0

    def __init__(self, setup=None):
        self._setup = setup
        self._use_setup = self.test_setup()

    def test_setup(self):
        """Determine whether we need a setup script, and which environment to use"""
        if self.setup and which(self.exe):
            # TODO: Add validation that the setup script works
            log.info('"%s" not detected, using supplied environment' % self.name)
            return True
        elif self.setup and which(self.exe) is None:
            log.info('%s detected, but using supplied environment instead.' % self.name + \
                     'Do not pass a "setup" argument to use local environment')
            return True
        elif self.setup is None and which(self.exe):
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

    @property
    def commandline(self):
        raise NotImplementedError("Subclasses should implement this!")

    @property
    def commandline_args(self):
        return shlex.split( self.commandline )

    def set_arguments(self, **kwargs):
        raise NotImplementedError("Subclasses should implement this!")

    def set_defaults(self):
        raise NotImplementedError("Subclasses should implement this!")

    def check_arguments(self):
        raise NotImplementedError("Subclasses should implement this!")

    def check_output(self):
        raise NotImplementedError("Subclasses should implement this!")

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

    def write_shell_script( self ):
        script_file = tempfile.NamedTemporaryFile(suffix='.sh', delete=False)
        with open( script_file, 'w') as handle:
            handle.write('source %s\n' % self.setup)
            handle.write( '%s\n' % self.commandline )
        return script_file

    def get_log_path( self, name ):
        log_name = '%s_%s.log' % (self._counter, name)
        return os.path.join( self._logs, log_name )

    def run( self ):
        log.info('Running %s' % self.name)
        if self.use_setup:
            log.info('Executing %s indirectly via Shell Script' % self.name)
        else:
            log.info('Executing %s directly via Subprocess' % self.name)
            print self.commandline
            output = self.run_subprocess()
            print output
        log.info('Finished Executing %s' % self.name)

    def run_subprocess_script(self):
        shell_script = self.write_shell_script()

    def run_subprocess(self):
        """Run a process directly via Subprocess and validate the output"""
        output = subprocess.check_output( self.commandline_args )
        return output

    def __call__(self, output_dir=None, **kwargs):
        self.counter += 1
        self.set_arguments( **kwargs )  # 1. Initialize
        self.set_defaults()             # 2. Set unspecified values to default
        self.check_arguments()          # 3. Sanity Check
        self.run()                      # 4. Run
        return self.check_output()      # 5. Validate