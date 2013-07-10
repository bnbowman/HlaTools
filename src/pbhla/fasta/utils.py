from string import maketrans

from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter

from pbhla.utils import make_rand_string, getbash, runbash

COMPLEMENT = maketrans('ACGT', 'TGCA')

def reverse_complement( fasta_record ):
    rev_seq = fasta_record.sequence[::-1]
    rev_com_seq = rev_seq.translate( COMPLEMENT )
    return FastaRecord( fasta_record.name,
                        rev_com_seq )

def write_fasta(fasta_records, output_file):
    with FastaWriter( output_file ) as handle:
        for record in fasta_records:
            handle.writeRecord( record )

def fasta_size(fasta):
    try:
        f = FastaReader(fasta)
        count=0
        for read in f:
            count+=1
        return count
    except:
	return None

def fasta_length( fasta ):
    try:
        f = FastaReader( fasta )
    except:
        return 0
    return max([len(read.sequence) for read in f])

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

def copy_fasta( fasta, destination, name=None ):
    with FastaWriter( destination ) as handle:
        for record in FastaReader( fasta ):
            if name:
                record._name = name
            handle.writeRecord( record )

def combine_fasta( fasta_files, destination ):
    with FastaWriter( destination ) as handle:
        for fasta in fasta_files:
            for record in FastaReader( fasta ):
                handle.writeRecord( record )
