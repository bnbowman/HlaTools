from collections import namedtuple

from pbcore.io.base import ReaderBase, getFileHandle
from pbhla.utils import BlasrM1, BlasrM4, BlasrM5

class BlasrReader( ReaderBase ):

    def __init__(self, f, filetype=None):
        self.file = getFileHandle(f, 'r')

        filetype = filetype or f.split('.')[-1]
        if filetype.lower() == 'm1':
            self.filetype = 'm1'
            self.datatype = BlasrM1
        elif filetype.lower() == 'm4':
            self.filetype = 'm4'
            self.datatype = BlasrM4
        elif filetype.lower() == 'm5':
            self.filetype = 'm5'
            self.datatype = BlasrM5
        else:
            raise TypeError("Invalid type to BlasrReader")

    def __iter__(self):
        try:
            for line in self.file:
                yield self.datatype._make( line.strip().split() )
        except:
            raise ValueError("Invalid Blasr entry of type %s" % self.filetype)
