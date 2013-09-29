#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

import os

from tempfile import NamedTemporaryFile
from pbhla.amplicon_analysis.sequence import (AmpliconAnalysisSequenceReader,
                                              AmpliconAnalysisSequenceWriter)
from pbhla.external.utils import align_best_reference
from pbhla.io.BlasrIO import BlasrReader



class ChimeraDetector(object):
    """
    A tool for detecting chimeric Amplicon Analysis consensus sequences
    """
    def __init__(self, query, mode="denovo", chunks=4):
        self._query = query
        self._sequences = self._parse_query_sequences()
        self._mode = mode
        self._chunks = chunks
        self._database = []

    def _parse_query_sequences(self):
        sequences = list(AmpliconAnalysisSequenceReader(self._query))
        return sorted(sequences, key=lambda s: s.num_reads, reverse=True)

    def run(self):
        """
        Search the input sequences for chimeric consensus reads
        """
        for record in self._sequences:
            # Cannot test for Chimeras without at least two references
            if len(self._database) < 2:
                self._database.append(record)
                continue
            # Write current database out to file
            db_file = _write_temp_fasta(self._database)

            # Chunk the query sequence and align the results to the DB
            chunks = list(_chunk_record(record, self._chunks))
            chunk_file = _write_temp_fasta(chunks)
            hits = _align_sequences(chunk_file, db_file)
            os.unlink(db_file)
            os.unlink(chunk_file)

            # Identify chimeric hits, adding non-chimeras to the database
            reference_hits = list(set([hit.tname for hit in hits]))
            if len(reference_hits) == 1:
                self._database.append(record)
            elif self._are_chimeric(hits):
                pass
            else:
                self._database.append(record)

def _chunk_record(record, chunks):
    """
    Slice a sequence record into multiple roughly-equal sized chunks
    """
    chunk_size = (len(record) / chunks) + 1
    for start in range(0, len(record), chunk_size):
        end = min(len(record), start+chunk_size)
        yield record[start:end]


def _write_temp_fasta(sequences):
    """
    Write the current sequence database out to file
    """
    temp = NamedTemporaryFile(suffix='.fasta', delete=False)
    with AmpliconAnalysisSequenceWriter(temp.name) as writer:
        for record in sequences:
            writer.write_fasta(record)
    return temp.name


def _align_sequences(query, reference):
    """
    Align one fasta file of sequences to another
    """
    temp = NamedTemporaryFile(suffix='.m1', delete=False)
    align_best_reference(query, reference, output=temp.name)
    hits = list(BlasrReader(temp.name))
    os.unlink(temp.name)
    return hits


if __name__ == "__main__":
    import sys

    query_file = sys.argv[1]

    ChimeraDetector(query_file).run()