from collections import namedtuple

from pbcore.io.base import ReaderBase, getFileHandle

BlasrM1 = namedtuple('BlasrM1', 'qname tname qstrand tstrand score pctsimilarity tstart tend tlength qstart qend qlength ncells')
BlasrM4 = namedtuple('BlasrM4', 'qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart tend tseqlength mapqv ncells clusterScore probscore numSigClusters')
BlasrM5 = namedtuple('BlasrM5', 'qname qlength qstart qend qstrand tname tlength tstart tend tstrand score nmat nmis nins ndel mapqv qstring astring tstring')

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

def add_header_to_m5( m5_file ):
    with open( m5_file ) as handle:
        lines = list( handle )
    with open( m5_file, 'w') as handle:
        handle.write('qname qlength qstart qend qstrand tname tlength tstart tend tstrand score nmat nmis nins ndel mapqv qstring astring tstring\n')
        for line in lines:
            handle.write( line )

def blasr_to_string( record ):
    if isinstance(record, BlasrM5):
        return "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (record.qname,
                                                                             record.qlength,
                                                                             record.qstart,
                                                                             record.qend,
                                                                             record.qstrand,
                                                                             record.tname,
                                                                             record.tlength,
                                                                             record.tstart,
                                                                             record.tend,
                                                                             record.tstrand,
                                                                             record.score,
                                                                             record.nmat,
                                                                             record.nmis,
                                                                             record.nins,
                                                                             record.ndel,
                                                                             record.mapqv,
                                                                             record.qstring,
                                                                             record.astring,
                                                                             record.tstring)
