#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

from pbcore.io.FastaIO import FastaRecord, FastaReader
from pbhla.io.utils import slice_2d

COLUMNS = 60

class FastaAlignment( object ):
    """
    A Class for representing a collection of FastaRecords as a single Multi-Sequence Alignment
    """

    def __init__(self, records):
        try:
            assert all( type(r) is FastaRecord for r in records )
            assert all( len(r.sequence) == len(records[0].sequence)
                        for r in records )
            self._records = records
        except AssertionError:
            raise ValueError("Invalid FASTA alignment data")

    @property
    def records(self):
        return self._records

    @property
    def size(self):
        return len(self.records[0].sequence)

    def __len__(self):
        return len(self.records)

    def __getitem__(self, args):
        # Return individual sequence Alignments if given Int
        rec_slice, seq_slice = slice_2d( args )
        records = self.records[rec_slice]
        sliced_records = [FastaRecord(r.name, r.sequence[seq_slice])
                             for r in records]
        filtered_records = [r for r in sliced_records
                               if len(set(r.sequence)) > 1]
        return FastaAlignment( filtered_records )

    def __iter__(self):
        for record in self.records:
            yield record

    def __str__(self):
        retval = ''
        for record in self.records:
            retval += ">%s\n" % record.name
            retval += "%s\n" % record.sequence
        return retval



def read_alignment( filename ):
    records = list(FastaReader( filename ))
    return FastaAlignment( records )

def write_fasta_alignment( filename, alignment ):
    with open( filename, 'w' ) as handle:
        for record in alignment:
            handle.write(">%s\n" % record.name)
            handle.write("%s\n" % record.sequence)

def write_stockholm_alignment( filename, alignment ):
    name_size = max([len(r.name) for r in alignment])
    with open( filename, 'w' ) as handle:
        handle.write("# STOCKHOLM 1.0\n\n")
        for start in xrange(0, alignment.size, COLUMNS):
            for record in alignment:
                handle.write("%s %s\n" % (record.name.ljust(name_size+1),
                                          record.sequence[start:start+COLUMNS]))
            handle.write("\n")
        handle.write('//')