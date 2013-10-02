from collections import namedtuple

from pbcore.io.base import ReaderBase, WriterBase, getFileHandle

blasr_m1_spec = 'qname tname qstrand tstrand score pctsimilarity tstart tend tlength qstart qend qlength ncells'
blasr_m4_spec = 'qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart tend tseqlength ' + \
                'mapqv ncells clusterScore probscore numSigClusters'
blasr_m5_spec = 'qname qlength qstart qend qstrand tname tlength tstart tend tstrand score nmat nmis nins ndel ' + \
                'mapqv qstring astring tstring'
BlasrM1 = namedtuple('BlasrM1', blasr_m1_spec)
BlasrM4 = namedtuple('BlasrM4', blasr_m4_spec)
BlasrM5 = namedtuple('BlasrM5', blasr_m5_spec)

class BlasrReader( ReaderBase ):

    def __init__(self, f, filetype=None):
        self._file = getFileHandle(f, 'r')

        filetype = filetype or f.split('.')[-1]
        if filetype.lower() == 'm1':
            self._filetype = 'm1'
            self._datatype = BlasrM1
        elif filetype.lower() == 'm4':
            self._filetype = 'm4'
            self._datatype = BlasrM4
        elif filetype.lower() == 'm5':
            self._filetype = 'm5'
            self._datatype = BlasrM5
        else:
            raise TypeError("Invalid type to BlasrReader")

    def __iter__(self):
        try:
            for line in self._file:
                entry = self._datatype._make(line.strip().split())
                if entry.qname == 'qname':
                    continue
                yield entry
        except:
            raise ValueError("Invalid Blasr entry of type %s" % self.filetype)



class BlasrWriter( WriterBase ):
    """
    A Class for writing out Blasr records
    """

    def writeHeader( self ):
        """
        Write a Blasr header out to the file handle 
        """
        self.file.write("qname tname qstrand tstrand score pctsimilarity tstart tend tlength qstart qend qlength ncells")
        self.file.write("\n")

    def write( self, record ):
        """
        Write a Blasr record out to the file handle
        """
        assert isinstance( record, BlasrM1 ) or isinstance( BlasrM5 )
        self.file.write( record_to_string( record ))
        self.file.write("\n")



def add_header_to_m5( m5_file ):
    with open( m5_file ) as handle:
        lines = list( handle )
    with open( m5_file, 'w') as handle:
        handle.write('qname qlength qstart qend qstrand tname tlength tstart tend tstrand score nmat nmis nins ndel mapqv qstring astring tstring\n')
        for line in lines:
            handle.write( line )

def record_to_string( record ):
    if isinstance(record, BlasrM1):
        return "%s %s %s %s %s %s %s %s %s %s %s %s %s" % (record.qname,
                                                           record.tname,
                                                           record.qstrand,
                                                           record.tstrand,
                                                           record.score,
                                                           record.pctsimilarity,
                                                           record.tstart,
                                                           record.tend,
                                                           record.tlength,
                                                           record.qstart,
                                                           record.qend,
                                                           record.qlength,
                                                           record.ncells)
    elif isinstance(record, BlasrM5):
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
