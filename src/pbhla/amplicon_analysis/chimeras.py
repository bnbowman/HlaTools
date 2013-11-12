#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'
__version__ = 0.3

import os
import logging

from tempfile import NamedTemporaryFile
from pbhla.amplicon_analysis.sequence import (AmpliconAnalysisSequenceReader,
                                              AmpliconAnalysisSequenceWriter)
from pbhla.external.utils import align_best_reference, multi_sequence_alignment
from pbhla.io.BlasrIO import BlasrReader
from pbhla.io.alignments import FastaAlignment
from pbhla.utils import valid_file

CHUNKS = 4
THRESHOLD = 0.28
BETA = 8
PSEUDOCOUNT = 1.4
MIN_COUNT = 20

log = logging.getLogger(__name__)

class ChimeraDetector(object):
    """
    A tool for detecting chimeric Amplicon Analysis consensus sequences
    """
    def __init__(self, debug=False,
                       prefix=None,
                       mode="denovo",
                       chunks=CHUNKS,
                       threshold=THRESHOLD,
                       beta=BETA,
                       pseudocount=PSEUDOCOUNT,
                       min_count=MIN_COUNT):
        log.info("Initializing ChimeraDetector v{0}".format(__version__))
        self._debug = debug

        self._mode = mode
        self._chunks = chunks
        self._threshold = threshold
        self._beta = beta
        self._pseudocount = pseudocount
        self._min_count = min_count
        self._database = {}
        self._chimeras = []

    @property
    def beta(self):
        return self._beta

    @property
    def threshold(self):
        return self._threshold

    @property
    def pseudocount(self):
        return self._pseudocount

    def _add_record(self, record):
        """
        Add a new sequence record to the
        """
        if record.name in self._database:
            msg = "Duplicate record (%s)" % record.name
            log.error( msg )
            raise KeyError( msg )
        log.info('Adding "%s" to the reference database\n' % record.name)
        self._database[record.name] = record

    def run(self, query, prefix=None):
        """
        Search the input sequences for chimeric consensus reads
        """
        if self._mode in ["d", "denovo"]:
            log.info("Running de novo Chimera detection tool")
            self._run_denovo_detection( query )
        elif self._mode in ["r", "reference"]:
            log.info('Running reference-based Chimera detection tool')
            pass
        else:
            msg = 'Mode must be either "denovo" or "reference"'
            log.error( msg )
            raise ValueError( msg )
        print [n for n in self._database]
        print [r.name for r in self._chimeras]
        self.write_output( query, prefix )

    def _run_denovo_detection(self, seq_file ):
        for record in sorted_sequences( seq_file ):
            log.info('Analyzing "%s"' % record.name)
            # Cannot test for Chimeras without at least two references
            if len(self._database) < 2:
                self._add_record( record )
                continue

            # Skip low-quality/coverage reads
            if record.num_reads < self._min_count:
                log.info("Skipping low-coverage sequence ({0} < {1})\n".format(record.num_reads,
                                                                               self._min_count))
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
            if hits is None:
                log.info("No hits found, probable new locus\n")
                self._add_record( record )
                continue

            # Identify chimeric hits, adding non-chimeras to the database
            reference_hits = list(set( [hit.tname for hit in hits] ))
            if len(reference_hits) <= 1:
                log.info('"%s" passed initial Chimera filter' % record.name)
                self._add_record( record )
                continue
            log.info('"%s" flagged as a possible Chimera with %s hits' % (record.name,
                                                                          len(reference_hits)))
            if len(reference_hits) == 2:
                references = [self._database[name] for name in reference_hits]
            elif len(reference_hits) > 2:
                best_hits = select_best_references(hits, 2)
                references = [self._database[name] for name in best_hits]

            if self.sequence_is_chimeric( record, references ):
                self._chimeras.append( record )
            else:
                self._add_record( record )

    def sequence_is_chimeric(self, query, references):
        """
        Determine whether a sequence in an Multi-Sequence alignment is Chimeric
        """
        # Create an MSA of the various sequences
        sequences = references + [query]
        alignment_file = _align_multiple_sequences(sequences)
        alignment = FastaAlignment( alignment_file )
        if self._debug:
            log.debug('Multi-sequence alignment for "%s" at %s' % (query.name,
                                                                   alignment_file))
        else:
            os.unlink( alignment_file )

        # Find which diffs align to which references
        votes = score_alignment_differences(query, references, alignment)
        #breakpoint, orientation = find_chimeric_breakpoint(votes)
        score = self.score_chimera(alignment, votes)
        if score > self.threshold:
            log.info("%s appears to be chimeric (%.2f)\n" % (query.name, score))
            return True
        else:
            log.info("%s appears to be non-chimeric (%.2f)" % (query.name, score))
            return False

    def score_chimera(self, alignment, votes):
        """
        Score chimeras at all possible break points
        """
        max_score = 0.0
        for breakpoint in sorted(votes):
            for left_ref, right_ref in [('A', 'B'), ('B', 'A')]:
                left = self.score_segment(votes, alignment.start, breakpoint, left_ref)
                right = self.score_segment(votes, breakpoint+1, alignment.end, right_ref)
                score = left*right
                max_score = max(score, max_score)
        return max_score

    def score_segment(self, votes, start, end, reference):
        yes = 0
        no = 0
        abstain = 0
        for position in votes:
            if position < start or position > end:
                continue
            if votes[position] == 'Abstain':
                abstain += 1
            elif votes[position] == reference:
                yes += 1
            else:
                no += 1
        return yes/((self.beta*(no+self.pseudocount))+abstain)

    def write_output(self, query, prefix):
        """
        Write the Chimeric and Non-Chimeric sequences out to separate files
        """
        prefix = prefix or '.'.join( query.split('.')[:-1] )
        suffix = query.split('.')[-1]
        bad_output = prefix + '.bad.' + suffix
        good_output = prefix + '.good.' + suffix

        log.info("Writing chimeric reads to {0}".format(bad_output))
        with AmpliconAnalysisSequenceWriter( bad_output ) as writer:
            for record in self._chimeras:
                writer.write_record( record )

        log.info("Writing non-chimeric reads to {0}".format(good_output))
        with AmpliconAnalysisSequenceWriter( good_output ) as writer:
            for n, record in self._database.iteritems():
                writer.write_record( record )


def score_alignment_differences(query, references, alignment):
    """
    Score diffs in alignment based on which reference they correspond to
    """
    A, B = references
    votes = {}
    for i in alignment.differences:
        if i < alignment.start or i > alignment.end:
            continue
        query_base = alignment[query.name][i]
        A_base = alignment[A.name][i]
        B_base = alignment[B.name][i]
        if query_base == A_base and query_base != B_base:
            votes[i] = 'A'
        elif query_base != A_base and query_base == B_base:
            votes[i] = 'B'
        else:
            votes[i] = 'Abstain'
    return votes

def select_best_references(hits, count):
    """
    Return the reference names for the top X Blasr hits by score
    """
    sorted_hits = sorted(hits, key=lambda x: float(x.pctsimilarity), reverse=True)
    first, rest = sorted_hits[0], sorted_hits[1:]
    second = [h for h in rest if h.tname != first.tname][0]
    return [first.tname, second.tname]

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
    if valid_file(temp.name):
        hits = list(BlasrReader(temp.name))
        os.unlink(temp.name)
        return hits
    os.unlink(temp.name)
    return None


def _align_multiple_sequences(sequences, output=None):
    """
    Align
    """
    temp_in = _write_temp_fasta( sequences )
    temp_out = NamedTemporaryFile(suffix=".afa", delete=False)
    multi_sequence_alignment( temp_in, temp_out.name )
    os.unlink( temp_in )
    return temp_out.name


def sorted_sequences( sequence_file ):
    """
    Yield successive reads in descending order of NumReads
    """
    seqs = list(AmpliconAnalysisSequenceReader( sequence_file ))
    for seq in sorted(seqs, key=lambda s: s.num_reads, reverse=True):
        yield seq

if __name__ == "__main__":
    import sys
    from pbhla.log import initialize_logger

    query_file = sys.argv[1]
    debug_flag = bool(sys.argv[2]) if len(sys.argv) > 2 else None

    initialize_logger( log_level='INFO' )
    cd = ChimeraDetector(debug=debug_flag)
    cd.run( query_file )
