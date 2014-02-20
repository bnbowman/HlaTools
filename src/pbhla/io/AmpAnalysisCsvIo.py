from collections import namedtuple

from pbcore.io.base import ReaderBase, WriterBase, getFileHandle

amp_analysis_definition = 'BarcodeName FastaName CoarseCluster Phase TotalCoverage Position Base QV Coverage HetScore'
AmpAnalysis = namedtuple('AmpAnalysis', amp_analysis_definition )

class AmpAnalysisReader( ReaderBase ):
    """
    A Class for parsing AmpliconAnalysis CSV files
    """

    def __iter__( self ):
        try:
            for line in self.file:
                entry = AmpAnalysis._make( line.strip().split(',') )
                if entry.BarcodeName == 'BarcodeName':
                    continue
                yield entry
        except:
            raise ValueError("Invalid AmpliconAnalysis entry")



class AmpAnalysisWriter( WriterBase ):
    """
    A Class for writing out AmpliconAnalysis records
    """

    def write_header( self ):
        """
        Write an AmpAnalysis header out to the file handle 
        """
        self.file.write( ','.join( amp_analysis_definition.split() ))
        self.file.write("\n")

    def write_record( self, record ):
        """
        Write an AmpAnalysis record out to the file handle
        """
        try:
            assert isinstance( record, AmpAnalysis )
            self.file.write( amp_analysis_to_string( record ))
            self.file.write("\n")
        except:
            raise ValueError("Invalid AmpliconAnalysis entry")



def amp_analysis_to_string( record ):
    """
    Convert an AmpAnalysis record to a string
    """
    assert isinstance( record, AmpAnalysis )
    data = [getattr(record, attr) for attr in amp_analysis_definition.split()]
    return ','.join( data )
