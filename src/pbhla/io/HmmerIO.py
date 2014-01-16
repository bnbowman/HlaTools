from pbcore.io.base import ReaderBase, getFileHandle

class HmmerDomainHit( object ):

    def __init__(self, *args):
        if len(args) == 1 and isinstance( args[0], str ):
            args = args[0].strip().split()
        try:
            assert len(args) == 23
            self.tname = args[0]
            self.qname = args[3]
            self.qlen = int(args[5])
            self.evalue = float(args[12])
            self.score = float(args[13])
            self.qstart = int(args[15])
            self.qend = int(args[16])
            self.tstart = int(args[17])
            self.tend = int(args[18])
        except:
            raise ValueError


class HmmerDomainReader( ReaderBase ):

    def __init__(self, f):
        self.file = getFileHandle(f, 'r')

    def __iter__(self):
        try:
            for line in self.file:
                if line.startswith('#'):
                    continue
                yield HmmerDomainHit( line )
        except:
            raise ValueError("Invalid Hmmer Domain Table entry")