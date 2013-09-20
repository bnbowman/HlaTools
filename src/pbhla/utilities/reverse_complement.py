import logging
from string import maketrans

from pbcore.io.FastaIO import FastaRecord
from pbcore.io.FastqIO import FastqRecord

COMPLEMENT = maketrans('ACGT', 'TGCA')

log = logging.getLogger(__name__)

def reverse_complement( record ):
    """
    Reverse complement a FastaRecord
    """
    rev_seq = record.sequence[::-1]
    rev_com_seq = rev_seq.translate( COMPLEMENT )
    if isinstance( record, FastaRecord ):
        return FastaRecord( record.name,
                            rev_com_seq )
    elif isinstance( record, FastqRecord ):
        rev_com_qual = record.qualityString[::-1]
        return FastqRecord( record.name,
                            rev_com_seq,
                            qualityString=rev_com_qual )
    else:
        msg = 'Record must be either Fasta or Fastq'
        log.error( msg )
        raise TypeError( msg )
