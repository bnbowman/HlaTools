import logging

from pbcore.io.FastaIO import FastaRecord, FastaReader, FastaWriter
from pbcore.io.FastqIO import FastqRecord, FastqReader, FastqWriter
from pbhla.utils import check_output_file, is_fasta, is_fastq

log = logging.getLogger(__name__)

def write_sequences( records, output_file ):
    """
    Write a sequence record, or list of records, out to file
    """
    if isinstance( records, list ) and all([isinstance( r, FastaRecord ) for r in records]):
        write_fasta( records, output_file )
    elif isinstance( records, list ) and all([isinstance( r, FastaRecord ) for r in records]):
        write_fastq( records, output_file )
    elif isinstance( records, FastaRecord ):
        write_fasta( [records], output_file )
    elif isinstance( records, FastqRecord ):
        write_fastq( [records], output_file )
    else:
        msg = "Input Record(s) type not recognized"
        log.error( msg )
        raise TypeError( msg )
    check_output_file( output_file )
    return output_file

def read_names( sequence_file ):
    # Open the sequence file with the appropriate reader
    if is_fasta( sequence_file ):
        reader = FastaReader( sequence_file )
    elif is_fastq( sequence_file ):
        reader = FastqReader( sequence_file )
    else:
        raise ValueError

    # Extract and return the sequence names
    return [r.name.strip().split()[0] for r in reader]

def write_fasta( records, output_file ):
    """
    Write a FastaRecord, or a list of FastaRecords, out to file
    """
    with FastaWriter( output_file ) as handle:
        for record in records:
            assert isinstance( record, FastaRecord )
            handle.writeRecord( record )
    check_output_file( output_file )
    return output_file

def write_fastq( records, output_file ):
    """
    Write a FastqRecord, or a list of FastqRecords, out to file
    """
    with FastqWriter( output_file ) as handle:
        for record in records:
            assert isinstance( record, FastqRecord )
            handle.writeRecord( record )
    check_output_file( output_file )
    return output_file

def combine_sequences( sequence_files, output_file ):
    """
    Combine a series of sequence files into one Fasta or Fastq
    """
    if all([is_fasta(f) for f in sequence_files]):
        combine_fasta( sequence_files, output_file )
    elif all([is_fastq(f) for f in sequence_files]):
        combine_fastq( sequence_files, output_file )
    else:
        msg = "Input files must be all Fasta or Fastq"
        log.error( msg )
        raise TypeError( msg )
    check_output_file( output_file )
    return output_file

def combine_fasta( sequence_files, output_file ):
    """
    Combine a series of sequence files into one Fasta
    """
    with FastaWriter( output_file ) as handle:
        for filename in sequence_files:
            try:
                for record in FastaReader( filename ):
                    handle.writeRecord( record )
            except:
                log.warn('Could not open "%s" as Fasta' % fasta)
    check_output_file( output_file )
    return output_file

def combine_fastq( sequence_files, output_file ):
    """
    Combine a series of sequence files into one Fastq
    """
    with FastqWriter( output_file ) as handle:
        for filename in sequence_files:
            try:
                for record in FastqReader( filename ):
                    handle.writeRecord( record )
            except:
                log.warn('Could not open "%s" as Fastq' % fasta)
    check_output_file( output_file )
    return output_file

def read_sequences( sequence_file ):
    """
    Parse a list of records from either a Fasta or Fastq file
    """
    if is_fasta( sequence_file ):
        return list( FastaReader( sequence_file ))
    elif is_fastq( sequence_file ):
        return list( FastqReader( sequence_file ))
    else:
        msg = 'Sequence file must be either Fasta or Fastq'
        log.error( msg )
        raise TypeError( msg )
