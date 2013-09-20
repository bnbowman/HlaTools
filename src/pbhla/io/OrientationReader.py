import logging

from pbhla.io.BlasrIO import BlasrReader

log = logging.getLogger(__name__)

class OrientationReader( object ):
    """
    A Class for parsing the orientation of reads from a Blasr alignment
    """

    def __init__( self, f ):
        if f.endswith('.m1') or f.endswith('.m4') or f.endswith('.m5'):
            self.orientations = _parse_orientation( f )
        else:
            msg = 'Invalid Orientation format (M1, M4, M5 valid)'
            log.error( msg )
            raise TypeError( msg )

    def is_forward( self, item ):
        """
        Return True if Forward, False if Reverse, otherwise None
        """
        try:
            if self.orientations[item] == 'forward':
                return True
            return False
        except:
            msg = 'Item not found! (%s)' % item
            log.error( msg )
            raise KeyError( msg )

    def is_reverse( self, item ):
        """
        Return False if Forward, True if Reverse, otherwise None
        """
        try:
            if self.orientations[item] == 'reverse':
                return True
            return False
        except:
            msg = 'Item not found! (%s)' % item
            log.error( msg )
            raise KeyError( msg )

    def __iter__( self ):
        return iter(self.orientations)

    def __contains__( self, item ):
        if item in self.orientations:
            return True
        return False



def _parse_orientation( filename ):
    """
    Parse the orientations of a list of sequences from a Blasr alignment file
    """
    orientations = {}
    for record in BlasrReader( filename ):
        if record.qname in orientations:
            msg = 'Duplicate record name! (%s)' % record.qname
            log.error( msg )
            raise ValueError( msg )
        if record.qstrand == record.tstrand:
            orientations[record.qname] = 'forward'
        else:
            orientations[record.qname] = 'reverse'
    return orientations
