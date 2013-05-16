from collections import namedtuple
 
entry = namedtuple('entry', 'qname flag rname pos mapq cigar rnext pnext tlen seq qual')

VALID_HD_TAGS = ['VN', 'SO']
VALID_SQ_TAGS = ['SN', 'LN', 'AS', 'M5', 'SP', 'UR']

REQUIRED_HD_TAGS = ['VN']
REQUIRED_SQ_TAGS = ['SN', 'LN']

class SamHeader( object ):
    def __init__(self, lines):
        self._version = None
        self._sort_order = None
        self._references = {}
        self._read_groups = []
        self._programs = []
        self._comments = []
        self._parse_input_lines( lines )

    def _parse_input_lines(self, lines):
        for line in lines:
            if line.startswith('@HD'):
                self._parse_header_line( line )
            elif line.startswith('@SQ'):
                self._parse_sequence_line( line )
            elif line.startswith('@RG'):
                self._parse_read_group_line( line )
            elif line.startswith('@PG'):
                self._parse_program_line( line )
            elif line.startswith('@CO'):
                self._parse_comment_line( line )
            else:
                msg = "Not a recognized header line: {0}".format( line )
                raise TypeError( msg )

    def _parse_header_line(self, line):
        if self._version:
            msg = "Only 1 header line allowed, but 2 detected"
            raise ValueError( msg )
        # Parse and validate the tags
        tags = tags_to_dictionary( line.strip().split()[1:] )
        validate_tags( tags, VALID_HD_TAGS, REQUIRED_HD_TAGS )
        # Set the appropriate variables
        self._version = tags['VN']
        if 'SO' in tags:
            self._sort_order = tags['SO']
        
    def _parse_sequence_line(self, line):
        tags = tags_to_dictionary( line.strip().split()[1:] )
        validate_tags( tags, VALID_SQ_TAGS, REQUIRED_SQ_TAGS )
        if tags['SN'] in self._references:
            msg = 'Sequence name "{0}" is duplicated!'.format(tags['SN'])
            raise ValueError( msg )
        tags['LN'] = int(tags['LN'])
        self._references[tags['SN']] = tags

    def _parse_read_group_line(self, line):
        pass

    def _parse_program_line(self, line):
        pass

    def _parse_comment_line(self, line):
        pass

    @property
    def version(self):
        return self._version

    @property
    def sort_order(self):
        return self._sort_order

    @property
    def references(self):
        return self._references

    @property
    def read_groups(self):
        return self._read_groups

    @property
    def program(self):
        return self._program

    @property
    def comments(self):
        return self._comments

class SamEntry( object ):
    def __init__(self, line):
        parts = line.strip().split()[:11]
        self.entry = entry._make(parts)
        self._pos = int(self.entry.pos)
        self._tlen = int(self.entry.tlen)

    @property
    def qname(self):
        return self.entry.qname

    @property
    def flag(self):
        return self.entry.flag

    @property
    def rname(self):
        return self.entry.rname

    @property
    def pos(self):
        return self._pos

    @property
    def mapq(self):
        return self.entry.mapq

    @property
    def cigar(self):
        return self.entry.cigar

    @property
    def rnext(self):
        return self.entry.rnext

    @property
    def pnext(self):
        return self.entry.pnext

    @property
    def tlen(self):
        return self._tlen

    @property
    def seq(self):
        return self.entry.seq

    @property
    def qual(self):
        return self.entry.qual

    @property
    def aend(self):
        return self.pos + self.tlen

class SamReader( object ):
    def __init__(self, f):
        self._file = open(f, "r")
        self._header = self.parse_header()
        self._file = open(f, "r") # Reset the file position

    def parse_header(self):
        header_lines = []
        line_start = 0
        for line in self._file:
            if line.startswith('@'):
                header_lines.append( line )
            else:
                break
        return SamHeader( header_lines )

    @property
    def header(self):
        return self._header

    @property
    def version(self):
        return self.header.version

    @property
    def sort_order(self):
        return self.header.sort_order

    @property
    def references(self):
        return self.header.references

    def close(self):
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        for line in self._file:
            if line.startswith('@'):
                continue
	    yield SamEntry(line)

#
# Utilities
#
def tags_to_dictionary( tags ):
    data_tags = {}
    for tag in tags:
        if tag[2] != ':':
            msg = 'Not a valid tag: "{0}"'.format(tag)
            raise TypeError( msg )
        tag_id, tag_value = tag[:2], tag[3:]
        data_tags[tag_id] = tag_value
    return data_tags

def validate_tags( tags, valid_tags, required_tags ):
    for tag in tags: # Check that all present tags are valid 
        if tag not in valid_tags:
            msg = 'Invalid tag "{0}" present'.format(tag)
            raise TypeError( msg )
    for tag in required_tags: # Check that all required tags are present
        if tag not in tags:
            msg = 'Required tag "{0}" not present'.format(tag)
            raise TypeError( msg )
