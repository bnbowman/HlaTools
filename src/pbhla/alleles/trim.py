#! /usr/bin/env python

import logging
from operator import itemgetter

from pbcore.io.FastqIO import FastqWriter, FastqRecord
from pbcore.io.FastaIO import FastaWriter, FastaRecord
from pbhla.fasta.utils import write_fasta
from pbhla.filenames import get_file_type
from pbhla.utils import check_output_file
from pbhla.sequences.utils import read_sequences

log = logging.getLogger()

def trim_alleles( input_file, output_file=None, trim=0 ):
    """Pick the top 2 Amplicon Analysis consensus seqs per group from a Fasta"""

    # If no trim or output file is specified, we can skip this module 
    if trim == 0 and output_file is None:
        log.info('No trimming necessary for "%s", skipping...' % input_file)
        return input_file

    # Set the output file if not specified
    output_file = output_file or _get_output_file( input_file )
    output_type = get_file_type( output_file )


    # Read the input sequences and trim each record 
    sequences = read_sequences( input_file )
    log.info('Trimming sequences by %s bp from each end' % trim)
    trimmed = _trim_sequences( sequences, trim )

    log.info('Writing the trimmed sequences out to %s' % output_file)
    _write_output( trimmed, output_file, output_type )
    return output_file

def _trim_sequences( records, trim ):
    """Trim X bases from each end of each sequence"""
    trimmed = []
    for record in records:
        if isinstance(record, FastaRecord):
            trimmed_record = FastaRecord( record.name, record.sequence[trim:-trim])            
        elif isinstance(record, FastqRecord):
            trimmed_record = FastqRecord( record.name, record.sequence[trim:-trim], record.quality[trim:-trim])
        else:
            raise TypeError("Only FastaRecord and FastqRecords support, not  '%s'" % type(record))
        trimmed.append( trimmed_record )
    return trimmed

def _write_output( records, output_file, output_type ):
    """Write the records out to file"""
    if output_type == 'fasta':
        write_fasta( records, output_file )
    else:
        with FastqWriter( output_file ) as writer:
            for record in records:
                writer.writeRecord( record )
        check_output_file( output_file )

def _get_output_file( input_file ):
    basename = '.'.join( input_file.split('.')[:-1] )
    file_type = get_file_type( input_file )
    return '%s.trimmed.%s' % (basename, file_type)

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    trim = int(sys.argv[3])
    
    trim_alleles( input_file, output_file, trim )
