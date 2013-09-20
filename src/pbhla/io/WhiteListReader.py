import logging

from pbcore.io.FastaIO import FastaReader
from pbcore.io.FastqIO import FastqReader

log = logging.getLogger(__name__)

class WhiteListReader( object ):
    """
    A Class for parsing reads to be White/Black-listed from
    """

    def __init__( self, f ):
        if f.endswith('.fa') or f.endswith('.fasta'):
            self.ids = _parse_fasta( f )
        elif f.endswith('.fq') or f.endswith('.fastq'):
            self.ids = _parse_fastq( f )
        elif f.endswith('.ids') or f.endswith('.txt'):
            self.ids = _parse_id_list( f )
        else:
            msg = 'Invalid WhiteList format (Fasta, Fastq, Ids and Txt valid)'
            log.error( msg )
            raise TypeError( msg )

    def __iter__( self ):
        return iter(self.ids)

    def __contains__( self, item ):
        if item in self.ids:
            return True
        return False



def _parse_fasta( filename ):
    return set([rec.name.split()[0] for rec in FastaReader( filename )])

def _parse_fastq( filename ):
    return set([rec.name.split()[0] for rec in FastqReader( filename )])

def _parse_id_list( filename ):
    return set([line.strip().split()[0] for line in open( filename ) if line])
