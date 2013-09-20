import re

HLA_PATTERN = re.compile('([A-Z]{1,3}\d?)[*_](\d{1,3})(:\d{1,3})?(:\d{1,3})?(:\d{1,3})?([NLSQ])?')

class HlaType( object ):
    """
    A Class for parsing and representing HLA types
    """

    def __init__( self, gene, 
                        field1, 
                        field2=None, 
                        field3=None, 
                        field4=None, 
                        suffix=None ):
        try:
            assert gene.isalnum()
            assert valid_field( field1 )
            self._gene = gene
            self._field1 = int(field1)

            if field2 is not None:
                assert valid_field( field2 )
                self._field2 = int(field2)
            else:
                self._field2 = None

            if field3 is not None:
                assert valid_field( field3 )
                self._field3 = int(field3)
            else:
                self._field3 = None

            if field4 is not None:
                assert valid_field( field4 )
                self._field4 = int(field4)
            else:
                self._field4 = None

            if suffix is not None:
                assert valid_suffix( suffix )
                self._suffix = suffix
            else:
                self._suffix = None
        except AssertionError:
            raise ValueError("Invalid HLA Type!")

    @property
    def gene( self ):
        return self._gene

    @property
    def field1( self ):
        return self._field1

    @property
    def field2( self ):
        return self._field2

    @property
    def field3( self ):
        return self._field3

    @property
    def field4( self ):
        return self._field4

    @property
    def suffix( self ):
        return self._suffix

    def __str__( self ):
        """
        Output a string representation of this HLA type
        """
        typing = 'HLA-%s*%02d' % (self.gene, self.field1)
        if self.field2 is not None:
            typing += ':%02d' % self.field2
        if self.field3 is not None:
            typing += ':%02d' % self.field3
        if self.field4 is not None:
            typing += ':%02d' % self.field4
        if self.suffix is not None:
            typing += self.suffix
        return typing

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return ( self.gene   == other.gene   and
                     self.field1 == other.field1 and 
                     self.field2 == other.field2 and
                     self.field3 == other.field3 and
                     self.field4 == other.field4 and
                     self.suffix == other.suffix )
        else:
            return False

    @classmethod
    def from_string( cls, s ):
        """
        Search a string for an HlaType and parse it
        """
        m = HLA_PATTERN.search( s )
        if m is None:
            return None
        field1 = None if m.group(2) is None else m.group(2)
        field2 = None if m.group(3) is None else m.group(3)[1:]
        field3 = None if m.group(4) is None else m.group(4)[1:]
        field4 = None if m.group(5) is None else m.group(5)[1:]
        return HlaType( gene=m.group(1),
                        field1=field1,
                        field2=field2,
                        field3=field3,
                        field4=field4,
                        suffix=m.group(6) )



def valid_field( field ):
    """
    Check whether a supplied variable is a valid HLA Type field
    """
    if isinstance( field, int ):
        if field > 0 and field < 1000:
            return True
    elif isinstance( field, str ):
        if field.isdigit() and len(field) <= 3:
            return True
    return False

def valid_suffix( suffix ):
    """
    Check whether a supplied variable is a valid HLA Type suffix
    """
    if isinstance( suffix, str ):
        if suffix in ['Q', 'S', 'L', 'N']:
            return True
    return False
