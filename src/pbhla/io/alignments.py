#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

from pbcore.io.FastaIO import FastaReader
from pbhla.io.BlasrIO import BlasrM5


class FastaAlignment(object):
    """
    Create an indexed alignment from an aligned Fasta file
    """
    NAME_SIZE = 18
    COLUMNS = 60

    def __init__(self, filename):
        try:
            self._records = list(FastaReader(filename))
            self._names = [r.name for r in self._records]
            self._sequences = [r.sequence for r in self._records]
            assert all([len(s) == len(self._sequences[0]) for s in self._sequences])
            self._differences = self._find_differences()
        except:
            raise ValueError("Invalid Fasta alignment data")

    def _find_differences(self):
        return [i for i in range(len(self.sequences[0])) if len(set(self[:,i])) != 1]

    @property
    def names(self):
        return self._names

    @property
    def sequences(self):
        return self._sequences

    @property
    def differences(self):
        return self._differences

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, item):
        if isinstance(item, int):
            if item < 0 or item > len(self):
                raise IndexError("The index (%d) is out of range" % item)
            return self.sequences[item]
        if isinstance(item, str):
            if item not in self.names:
                raise KeyError("")
            return self.sequences[self.names.index(item)]
        elif isinstance(item, slice):
            start, stop, step = item.indices(len(self))
            return self.sequences[start:stop:step]
        elif isinstance(item, tuple):
            assert len(item) == 2, "Alignment only has two dimensions"
            x, y = item
            return [sequence[y] for sequence in self[x]]
        else:
            raise TypeError("Index must be Int, Slice or Tuple")

    def __str__(self):
        retval = ''
        for i in range(0, len(self), self.COLUMNS):
            for n, s in zip(self.names, self.sequences):
                retval += '%s %s\n' % (n, s[i:i+self.COLUMNS])
            retval += '\n'
        return retval

class BlasrAlignment(object):
    """
    Create an indexed alignment from a Blasr M5 record
    """
    NAME_SIZE = 12
    COLUMNS = 60

    def __init__(self, record):
        try:
            assert isinstance(record, BlasrM5)
            assert len(record.qstring) == len(record.astring)
            assert len(record.qstring) == len(record.astring)
            self._query_name = trim_string(record.qname, self.NAME_SIZE)
            self._template_name = trim_string(record.tname, self.NAME_SIZE)
            self._query = record.qstring
            self._alignment = record.astring
            self._template = record.tstring
        except:
            raise ValueError("Invalid Blasr alignment data")

    @property
    def query_name(self):
        return self._query_name

    @property
    def template_name(self):
        return self._template_name

    @property
    def query(self):
        return self._query

    @property
    def alignment(self):
        return self._alignment

    @property
    def template(self):
        return self._template

    @property
    def differences(self):
        return [i for i in range(len(self)) 
                        if self.alignment[i] == '*']

    def __len__(self):
        return len(self.query)

    def __str__(self):
        retval = ''
        for i in range(0, len(self), self.COLUMNS):
            retval += '%s %s\n' % (self.query_name.ljust(self.NAME_SIZE),
                                   self.query[i:i+self.COLUMNS])
            retval += '%s %s\n' % (''.ljust(self.NAME_SIZE),
                                   self.alignment[i:i+self.COLUMNS])
            retval += '%s %s\n\n' % (self.template_name.ljust(self.NAME_SIZE),
                                     self.template[i:i+self.COLUMNS])
        return retval


def trim_string(s, size):
    if len(s) >= size:
        return s[:size]
    return s.ljust(size)

if __name__ == '__main__':
    import sys

    a = FastaAlignment(sys.argv[1])