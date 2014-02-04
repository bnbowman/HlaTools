#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

from pbhla.external.ExternalTool import ExternalTool
from pbhla.utils import validate_file

#Test

class ToAmosTool( ExternalTool ):
    """A simply wrapper for the toAmos data conversion tool"""

    def __init__(self, setup=None):
        super(ToAmosTool, self).__init__( setup=setup )
        self.output = None
        self.fasta_arg = None
        self.fastq_arg = None
        self.output_arg = None

    @property
    def name(self):
        return 'toAmos'

    @property
    def exe(self):
        return self.name

    @property
    def commandline(self):
        """Convert the supplied arguments into a command-line call"""
        cmd = "{exe} {s} {q} {o}"
        cmd = cmd.format(exe=self.exe, s=self.fasta_arg, q=self.fastq_arg, o=self.output_arg)
        return cmd

    def set_arguments(self, **kwargs):
        for key, value in kwargs.iteritems():
            if key in ['s', 'fasta']:
                self.fasta_arg = '-s %s' % value
            elif key in ['q', 'fastq']:
                self.fastq_arg = '-q %s' % value
            elif key in ['o', 'output']:
                self.output = value
                self.output_arg = '-o %s' % value
            else:
                raise ValueError('"%s" is not a valid argument' % key)

    def set_defaults(self):
        if not self.output_arg:
            if self.fasta_arg:
                self.output_arg = '-o %s' % 'blah'
            elif self.fastq_arg:
                pass
            else:
                raise ValueError

    def check_arguments(self):
        assert self.fasta_arg or self.fastq_arg
        assert self.output_arg and self.output

    def check_output(self):
        return validate_file( self.output )