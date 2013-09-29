import logging, tempfile
from string import maketrans

from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter
from pbcore.io.FastqIO import FastqRecord
from pbhla.utils import check_output_file

COMPLEMENT = maketrans('ACGT', 'TGCA')
log = logging.getLogger(__name__)

def write_temp_fasta( record ):
    """
    Write a temporary Fasta file
    """
    temp = tempfile.NamedTemporaryFile( suffix='.fasta', delete=False )
    if isinstance( record, FastaRecord ):
        write_fasta( record, temp.name )
    elif isinstance( record, FastqRecord ):
        fasta = FastaRecord( record.name, record.sequence )
        write_fasta( fasta, temp.name )
    else:
        msg = 'Sequence record must be either Fasta or Fastq'
        log.error( msg )
        raise TypeError( msg )
    return temp

def reverse_complement( fasta_record ):
    """
    Reverse complement a FastaRecord
    """
    rev_seq = fasta_record.sequence[::-1]
    rev_com_seq = rev_seq.translate( COMPLEMENT )
    return FastaRecord( fasta_record.name,
                        rev_com_seq )

def write_fasta( fasta_records, output_file):
    """
    Write a FastaRecord, or list of records, out to file
    """
    with FastaWriter( output_file ) as handle:
        if isinstance( fasta_records, FastaRecord ):
            handle.writeRecord( fasta_records )
        elif isinstance( fasta_records, list):
            for record in fasta_records:
                handle.writeRecord( record )
        else:
            msg = "Input Record(s) type not recognized"
            log.error( msg )
            raise TypeError( msg )
    check_output_file( output_file )

def fasta_size(fasta):
    """
    Count the number of sequences in a Fasta
    """
    try:
        f = FastaReader(fasta)
        count=0
        for read in f:
            count+=1
        return count
    except:
	return None

def fasta_length( fasta ):
    """
    Return the maximum sequence length in a Fasta file
    """
    try:
        f = FastaReader( fasta )
    except:
        return 0
    return max([len(read.sequence) for read in f])

def extract_names( fasta ):
    """
    Extract all of the names from a Fasta file
    """
    names = []
    for record in FastaReader( fasta ):
        name = record.name.split()[0]
        names.append( name )
    return names

def extract_sequence(fasta, names):
    f = FastaReader(fasta)
    if isinstance(names, str):
        for r in f:
            if r.name == names:
                return r.sequence
    elif isinstance(names, list):
        output=[]
        for r in f:
            if r.name in names:
                output.append(r)
        return output

def is_fasta( filename ):
    if filename.endswith('.fa') or filename.endswith('.fasta'):
        return True
    return False

def copy_fasta( fasta, destination, name=None ):
    with FastaWriter( destination ) as handle:
        for record in FastaReader( fasta ):
            if name:
                record._name = name
            handle.writeRecord( record )

def combine_fasta( fasta_files, destination ):
    with FastaWriter( destination ) as handle:
        for fasta in fasta_files:
            try:
                for record in FastaReader( fasta ):
                    handle.writeRecord( record )
            except:
                log.warn('Could not open "%s" as Fasta' % fasta)
    check_output_file( destination )
