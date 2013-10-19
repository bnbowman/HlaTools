#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

from pbcore.io.base import ReaderBase, WriterBase

class FofnReader( ReaderBase ):
    """
    A Class for reading file names from a FOFN
    """

    def __iter__(self):
        for line in self.file:
            filename = line.strip()
            if not filename:
                continue
            try:
                assert len(filename.split()) == 1
            except AssertionError:
                raise ValueError("Invalid FOFN entry (%s)" % filename)
            yield filename


class FofnWriter( WriterBase ):
    """
    A Class for writing file names to a FOFN
    """

    def write( self, filename ):
        """
        Write a FOFN record out to the file handle
        """
        assert isinstance( filename, str )
        assert len(filename.split()) == 1
        self.file.write( filename + "\n" )