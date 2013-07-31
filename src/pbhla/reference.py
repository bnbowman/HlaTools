import os, re, csv, logging

from pbhla.fasta.utils import write_fasta

log = logging.getLogger()

def create_reference_fasta( fofn_file, output_file ):
    log.info("Creating a reference Fasta from a FOFN of sequence files")
    log.debug("\tInput File:\t%s" % fofn_file)
    log.debug("\tOutput File:\t%s" % output_file)

    records = _parse_sequence_records( fofn_file )
    log.info("Found %s reference sequence records" % len(records))
    write_fasta( records, output_file )
    log.info("Finished creating reference Fasta file")

def _parse_sequence_records( fofn_file ):
    records = []
    with open( fofn_file, 'r') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            filename, locus = line.strip().split()
            records += list( FastaReader( filename ) )
    return records

def read_reference_metadata( fofn_file ):
    metadata = {}
    with open( fofn_file, 'r') as handle:
        for line in handle:
            if not line.startswith('#'):
                continue
            key, value = line[1:].strip().split()
            metadata[key] = value
    return metadata
