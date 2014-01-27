#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

from pbhla.external.ExternalTool import ExternalTool
from pbhla.utils import validate_file

class Minimus2Tool( ExternalTool ):
    """A simply wrapper for the toAmos data conversion tool"""

    def __init__(self, setup=None):
        super(Minimus2Tool, self).__init__( setup=setup )
        self.output = None
        self.prefix_arg = None

    @property
    def name(self):
        return 'minimus2'

    @property
    def exe(self):
        return self.name

    @property
    def commandline(self):
        """Convert the supplied arguments into a command-line call"""
        cmd = "{exe} {p}"
        cmd = cmd.format(exe=self.exe, p=self.prefix_arg)
        return cmd

    def set_arguments(self, **kwargs):
        for key, value in kwargs.iteritems():
            if key in ['a', 'afg']:
                self.prefix_arg = '.'.join( value.split('.')[:-1] )
            elif key in ['p', 'prefix']:
                self.prefix_arg = value
            else:
                raise ValueError('"%s" is not a valid argument' % key)

    def set_defaults(self):
        self.output = self.prefix_arg + '.fasta'

    def check_arguments(self):
        assert self.prefix_arg
        assert not validate_file(self.output)

    def check_output(self):
        return validate_file( self.output )