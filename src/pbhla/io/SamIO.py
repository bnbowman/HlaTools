from collections import namedtuple
 
entry = namedtuple('entry', 'qname flag rname pos mapq cigar rnext pnext tlen seq qual')

class SamReader:
    def __init__(self, f):
        self.file = open(f, "r")

    def close(self):
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        for line in self.file:
            if len(line) == 0: 
                continue
            if line.startswith('@'): 
                continue
            line = line.strip().split()[:11]
	    yield entry._make(line)
