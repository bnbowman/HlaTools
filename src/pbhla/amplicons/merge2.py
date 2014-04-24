#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import logging
import tempfile

from pbcore.io import FastaWriter, FastqReader, FastaRecord, FastqRecord

from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import fasta_size, read_fasta_dict, write_fasta
from pbhla.utils import read_fastq_dict, consensus_filetype, write_fastq
from pbhla.io.BlasrIO import BlasrReader

log = logging.getLogger()
log.setLevel( logging.INFO )

def merge_amplicons( sequence_5p, sequence_3p, reference, output ):
    file_list = [sequence_5p, sequence_3p]
    filetype = consensus_filetype( file_list )
    sequences = read_sequences( filetype, file_list )

    alignment_left  = align_amplicons( sequence_5p, reference )
    alignment_right = align_amplicons( sequence_3p, reference )
    positions_left  = parse_alignment_positions( alignment_left )
    positions_right = parse_alignment_positions( alignment_right )

    pairs = pair_sequences( positions_left, positions_right )
    check_overlap( pairs )

    merged = merge_sequences( filetype, sequences, pairs )
    write_sequences( filetype, merged, output )

def check_overlap( pairs ):
    for left, right in pairs:
        if right['tstart'] > left['tend']:
            raise ValueError('Pair from "%s" do not overlap' % left['ref'])
    return True

def get_major_allele( reference ):
    locus_string = reference.split('_')[1]
    major_allele = locus_string.split(':')[0]
    return major_allele

def pair_sequences( positions_left, positions_right ):
    left_refs  = sorted([p['ref'] for p in positions_left])
    right_refs = sorted([p['ref'] for p in positions_right])

    if len(left_refs) != len(right_refs):
        raise ValueError("Reference Mismatch")

    if left_refs == right_refs:
        pairs = []
        for ref in left_refs:
            left =  [p for p in positions_left  if p['ref'] == ref][0]
            right = [p for p in positions_right if p['ref'] == ref][0]
            pairs.append( (left, right) )
        return pairs

    left_majors = sorted([get_major_allele(r) for r in left_refs])
    right_majors = sorted([get_major_allele(r) for r in right_refs])

    if left_majors == right_majors:
        pairs = []
        for ref in left_majors:
            left =  [p for p in positions_left  if get_major_allele(p['ref']) == ref][0]
            right = [p for p in positions_right if get_major_allele(p['ref']) == ref][0]
            pairs.append( (left, right) )
        return pairs

    raise ValueError("Reference Mismatch")


def write_sequences( filetype, merged, output ):
    if filetype == 'fasta':
        write_fasta( merged, output )
    elif filetype == 'fastq':
        write_fastq( merged, output )
    else:
        raise ValueError

def merge_sequences( filetype, sequences, pairs ):
    if filetype == 'fastq':
        return merge_fastq_sequences( sequences, pairs )
    elif filetype == 'fasta':
        return merge_fasta_sequences( sequences, pairs )
    else:
        raise ValueError

def merge_fasta_sequences( sequences, pairs ):
    merged = []
    for i, part_list in enumerate( positions ):
        name = "Merged%s_NumReads100" % (i+1)
        sequence = ''
        for part in part_list:
            source = part['name']
            start = part['start']
            end = part['end']
            sequence += sequences[source].sequence[start:end]
        record = FastaRecord(name=name, sequence=sequence)
        merged.append(record)
    return merged

def merge_fastq_sequences( sequences, pairs ):
    merged = []
    for i, pair in enumerate( pairs ):
        left, right = pair
        left_rec = sequences[left['name']]
        right_rec = sequences[right['name']]
        name = "Merged%s_NumReads100" % (i+1)
        overlap = left['tend'] - right['tstart']
        qoverlap = left['qstring'][-overlap:]
        toverlap = left['tstring'][-overlap:]
        indel_correction = toverlap.count('-') - qoverlap.count('-')
        end = left['qend'] - overlap - indel_correction
        sequence = left_rec.sequence[:end] + right_rec.sequence
        quality = left_rec.qualityString[:end] + right_rec.qualityString
        record = FastqRecord(name=name, sequence=sequence, qualityString=quality)
        merged.append( record )
    return merged

def read_sequences( filetype, file_list ):
    if filetype == 'fastq':
        return read_fastq_dict( file_list )
    elif filetype == 'fasta':
        return read_fasta_dict( file_list )
    else:
        raise ValueError

def parse_alignment_positions( alignment_file ):
    positions = []
    for hit in BlasrReader( alignment_file ):
        position = {'name':hit.qname,
                    'ref': hit.tname,
                    'qstart':int(hit.qstart),
                    'qend':int(hit.qend),
                    'tstart':int(hit.tstart),
                    'tend':int(hit.tend),
                    'qstring':hit.qstring,
                    'tstring':hit.tstring}
        positions.append( position )
    return positions

def write_temp_fasta( fastq_file ):
    """
    Write a temporary Fasta file from a Fastq
    """
    temp = tempfile.NamedTemporaryFile( suffix='.fasta', delete=False )
    with FastaWriter( temp.name ) as handle:
        for record in FastqReader( fastq_file ):
            temp_record = FastaRecord( record.name, record.sequence )
            handle.writeRecord( temp_record )
    return temp

def align_amplicons( sequence_file, reference ):
    sequence_root = '.'.join( sequence_file.split('.')[:-1] )
    output_file = sequence_root + '.m5'
    blasr_args = {'bestn': 1,
                  'out': output_file,
                  'm': 5,
                  'noSplitSubreads': True}
    return run_blasr( sequence_file, reference, blasr_args, verbose=True )

if __name__ == '__main__':
    import sys

    sequence_5p = sys.argv[1]
    sequence_3p = sys.argv[2]
    reference = sys.argv[3]
    output_file = sys.argv[4]

    logging.basicConfig( stream=sys.stdout )
    merge_amplicons( sequence_5p, sequence_3p, reference, output_file )