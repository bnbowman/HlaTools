#! /usr/bin/env python

import logging
from tempfile import NamedTemporaryFile

from pbcore.io.FastaIO import FastaReader, FastaRecord
from pbcore.io.FastqIO import FastqReader, FastqRecord
from pbhla.filenames import get_file_type
from pbhla.io.BlasrIO import BlasrReader
from pbhla.io.utils import parse_locus_dict, get_output_file, get_temp_fasta, write_records
from pbhla.io.HmmerIO import HmmerDomainReader
from pbhla.external.commandline_tools import run_hmmsearch
from pbhla.external.utils import get_alignment_file

log = logging.getLogger()

def cdna_from_file( input_file, hmm_fofn, output=None,
                                          reference=None,
                                          alignment=None ):
    """
    Extract the cDNA sequences from a mixed Fasta or Fastq
    """
    # Check the input files, and align the input file if needed
    alignment_file = get_alignment_file( input_file, reference, alignment )
    output_file = output or get_output_file( input_file, 'cDNA' )

    # Prepare the Fasta by orienting and subsetting it
    records = _parse_input_records( input_file )
    hmms = parse_locus_dict( hmm_fofn )
    loci = _parse_loci( alignment_file )

    # Compose and output the records
    cdna_records = list( cdna_from_records( records, loci, hmms ))
    write_records( cdna_records, output_file )
    return output_file

def cdna_from_records( records, loci, hmms ):
    """
    Extract the cDNA sequences from a list of Records
    """
    for record in records:

        # First identify the locus of a given record
        name = record.name.split()[0]
        try:
            locus = loci[name]
        except KeyError:
            msg = 'No HLA locus associated with "%s" - skipping' % name
            log.warn( msg )
            #raise UserWarning( msg )
            continue

        # Find the HMM model appropriate for the locus
        try:
            hmm = hmms[locus]
        except KeyError:
            msg = 'No exonic HMM for %s - skipping' % locus
            log.warn( msg )
            #raise UserWarning( msg )
            continue

        # Extract the exons and make the cDNA
        yield cdna_from_record( record, hmm )

def cdna_from_record( record, hmm ):
    """
    Extract the cDNA sequence from a single record
    """
    # Write the record to file and search
    temp_fasta = get_temp_fasta( record )
    temp_hmm_output = NamedTemporaryFile()
    hmmsearch_args = {'domE': 1e-10,
                      'domtblout': temp_hmm_output.name}
    run_hmmsearch( hmm, temp_fasta.name, hmmsearch_args )
    temp_fasta.close()

    #with open( temp_hmm_output.name ) as handle:
    #    for line in handle:
    #        print line.strip()
    # Parse the exons from the HMM output
    exons = parse_exon_slices( temp_hmm_output.name )
    temp_hmm_output.close()

    # Extract and return the cDNA record
    cdna_record = _multislice_record( record, exons )
    return cdna_record

def _multislice_record(record, slices):
    """
    Slice and combine multiple regions from a Fasta or Fastq
    """
    sliced_records = [_slice_record(record, s) for name, s in slices]
    sequence = ''.join([r.sequence for r in sliced_records])
    if isinstance(record, FastaRecord):
        return FastaRecord(record.name, sequence)
    elif isinstance(record, FastqRecord):
        quality_str = ''.join([r.qualityString for r in sliced_records])
        return FastqRecord(record.name, sequence, qualityString=quality_str)
    else:
        msg = 'Invalid sequence record type'
        log.error( msg )
        raise TypeError( msg )

def _slice_record(record, slice):
    """
    Slice a region out of a Fasta or Fastq record
    """
    sequence = record.sequence[slice]
    if isinstance(record, FastaRecord):
        return FastaRecord(record.name, sequence)
    elif isinstance(record, FastqRecord):
        quality = record.quality[slice]
        return FastqRecord(record.name, sequence, quality)
    else:
        msg = 'Invalid sequence record type'
        log.error( msg )
        raise TypeError( msg )

def parse_exon_slices( hmm_output_file ):
    starts = {}
    ends = {}
    exons = []
    for hit in HmmerDomainReader( hmm_output_file ):
        start = hit.tstart - hit.qstart
        end = hit.tend + (hit.qlen - hit.qend)
        exons.append( (hit.qname, slice(start, end)) )
        #starts[hit.qname] = min( start, starts.get(hit.qname, 10000))
        #ends[hit.qname] = max( end, ends.get(hit.qname, 0))
    #exons = [(name, slice(starts[name], ends[name])) for name in starts]
    return sorted(exons, key=lambda x: x[0])

def _parse_input_records( input_file ):
    """
    Parse the input sequence records with the appropriate pbcore Reader
    """
    input_type = get_file_type( input_file )
    if input_type == 'fasta':
        return list( FastaReader( input_file ))
    elif input_type == 'fastq':
        return list( FastqReader( input_file ))
    else:
        msg = 'Input file must be either Fasta or Fastq'
        log.error( msg )
        raise ValueError( msg )

def _parse_loci( blasr_file ):
    """
    Parse the likely locus of sequences from a Blasr file
    """
    locus_calls = {}
    for entry in BlasrReader( blasr_file ):
        if entry.tname == 'tname':
            continue

        # Parse the locus from either Tokai or IMGT references
        reference = entry.tname.split('*')[0]
        if reference.startswith('HLA-'):
            locus = reference[-1]
        else:
            locus = reference.split('_')[1]

        # Save the Locus/Sequence pair unless duplicate
        if entry.qname in locus_calls:
            msg = 'Duplicate sequence name found "%s"!' % entry.qname
            log.error( msg )
            raise ValueError( msg )
        else:
            locus_calls[entry.qname] = locus
    return locus_calls

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO, stream=sys.stdout )

    input_file = sys.argv[1]
    hmm_fofn = sys.argv[2]
    reference = sys.argv[3]

    cdna_from_file( input_file, hmm_fofn,
                                reference=reference )