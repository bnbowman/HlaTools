__author__ = 'bbowman@pacificbiosciences.com'

import os
import logging

from tempfile import NamedTemporaryFile
from pbhla.amplicon_analysis.sequence import (AmpliconAnalysisSequenceReader,
                                              AmpliconAnalysisSequenceWriter)
from pbhla.external.utils import align_best_reference, multi_sequence_alignment
from pbhla.io.BlasrIO import BlasrReader
from pbhla.io.alignments import FastaAlignment


log = logging.getLogger(__name__)

class ChimeraDetector(object):
    """
    A tool for detecting chimeric Amplicon Analysis consensus sequences
    """
    def __init__(self, query, mode="denovo", chunks=4):
        log.info("Initializing ChimeraDetector")
        self._query = query
        self._sequences = self._parse_query_sequences()
        self._mode = mode
        self._chunks = chunks
        self._database = {}
        self._chimeras = []

    def _parse_query_sequences(self):
        """
        Parse and sort the query sequences by NumReads
        """
        sequences = list(AmpliconAnalysisSequenceReader(self._query))
        return sorted(sequences, key=lambda s: s.num_reads, reverse=True)

    def _add_record(self, record):
        """
        Add a new sequence record to the
        """
        if record.name in self._database:
            msg = "Duplicate record (%s)" % record.name
            log.error( msg )
            raise KeyError( msg )
        log.info('Adding "%s" to the reference database' % record.name)
        self._database[record.name] = record

    def run(self):
        """
        Search the input sequences for chimeric consensus reads
        """
        if self._mode in ["d", "denovo"]:
            log.info("Running de novo Chimera detection tool")
            self._run_denovo_detection()
        elif self._mode in ["r", "reference"]:
            log.info('Running reference-based Chimera detection tool')
            pass
        else:
            msg = 'Mode must be either "denovo" or "reference"'
            log.error( msg )
            raise ValueError( msg )
        print [n for n in self._database]
        print [r.name for r in self._chimeras]

    def _run_denovo_detection(self):
        for record in self._sequences:
            log.info('Analyzing "%s"' % record.name)
            # Cannot test for Chimeras without at least two references
            if len(self._database) < 2:
                self._add_record( record )
                continue

            # Write current database out to file
            db_records = list( self._database.itervalues() )
            db_file = _write_temp_fasta( db_records )

            # Chunk the query sequence and align the results to the DB
            chunks = list( _chunk_record( record, self._chunks ))
            chunk_file = _write_temp_fasta( chunks )
            hits = _align_sequences( chunk_file, db_file )
            os.unlink( db_file )
            os.unlink( chunk_file )

            # Identify chimeric hits, adding non-chimeras to the database
            reference_hits = list(set( [hit.tname for hit in hits] ))
            if len(reference_hits) <= 1:
                log.info('"%s" passed initial Chimera filter' % record.name)
                self._add_record( record )
                continue
            elif len(reference_hits) == 2:
                references = [self._database[name] for name in reference_hits]
            elif len(reference_hits) > 2:
                best_hits = select_best_references(hits, 2)
                references = [self._database[name] for name in best_hits]

            log.info('"%s" flagged as a possible Chimera' % record.name)
            if sequence_is_chimeric( record, references ):
                log.info('"%s" appears Chimeric, excluding...' % record.name)
                self._chimeras.append( record )
            else:
                log.info('"%s" appears non-Chimeric' % record.name)
                self._add_record( record )

def sequence_is_chimeric(query, references):
    """
    Determine whether a sequence in an Multi-Sequence alignment is Chimeric
    """
    sequences = references + [query]
    alignment_file = _align_multiple_sequences(sequences)
    alignment = FastaAlignment( alignment_file )
    print alignment
    print str(alignment)
    print alignment.differences
    return True

def select_best_references(hits, count):
    """
    Return the reference names for the top X Blasr hits by score
    """
    sorted_hits = sorted(hits, lambda x: x.score)
    selected_hits = sorted_hits[:count]
    return [hit.tname for hit in selected_hits]

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


def _align_multiple_sequences(sequences, output=None):
    """
    Align
    """
    temp_in = _write_temp_fasta( sequences )
    temp_out = NamedTemporaryFile(suffix=".afa", delete=False)
    multi_sequence_alignment( temp_in, temp_out.name )
    os.unlink( temp_in )
    return temp_out.name


if __name__ == "__main__":
    import sys
    logging.basicConfig( level=logging.INFO )

    query_file = sys.argv[1]

    ChimeraDetector(query_file).run()
