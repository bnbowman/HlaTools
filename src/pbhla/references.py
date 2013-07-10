import os, re, csv, logging

from pbcore.io.FastaIO import FastaReader
from pbhla.io.BlasrIO import BlasrReader, blasr_to_string
from pbhla.io.SamIO import SamReader
from pbhla.fasta.utils import write_fasta

log = logging.getLogger()

def create_fofn_reference( fofn_file ):
    log.info('Reading Locus References from "{0}"'.format(fofn_file))
    results = {}
    with open(fofn_file, 'r') as handle:
        for line in handle:
            fasta_path, locus = line.strip().split()
            fasta_file = os.path.basename( fasta_path )
            log.info('Reading "{0}" sequences from "{1}"'.format(locus, fasta_file))
            for record in FastaReader( fasta_path ):
                name = record.name.split()[0]
                if not re.search('__', name):
                    name = name.split('_')[0]
                results[name] = locus
    log.info('Finished reading Locus References')
    return results

def create_m1_reference( m1_file, reference=None ):
    log.info('Parsing Blasr M1 results from "{0}"'.format( m1_file ))
    results = {}
    for record in BlasrReader( m1_file ):
        if record.qname in results:
            msg = 'Duplicate sequence ids found! "{0}"'.format( record.qname )
            log.info( msg )
            raise KeyError( msg )
        if reference:
            results[record.qname] = reference[record.tname]
        else:
            results[record.qname] = record.tname
    log.info('Finished reading Blasr results')
    return results

def create_m5_reference( m5_file ):
    log.info('Parsing Blasr M5 results from "{0}"'.format( m5_file ))
    results = {}
    diffs = {}
    for record in BlasrReader( m5_file ):
        diff_count = int(record.nmis) + int(record.nins) + int(record.ndel)
        if record.qname not in diffs:
            results[record.qname] = record.tname
            diffs[record.qname] = diff_count
        elif diffs[record.qname] > diff_count:
            results[record.qname] = record.tname
            diffs[record.qname] = diff_count
    log.info('Finished reading Blasr results')
    return results

def create_sam_reference( sam_file, reference=None ):
    log.info('Parsing SAM alignments from "{0}"'.format(sam_file))
    results = {}
    for record in SamReader(sam_file):
        if record.qname in results:
            msg = 'Duplicate sequence ids found! "{0}"'.format( record.qname )
            log.info( msg )
            raise KeyError( msg )
        if reference:
            results[record.qname] = reference[record.rname]
        else:
            results[record.qname] = record.rname
    log.info('Finished reading SAM file results')
    return results

def filter_m5_file( m5_file, filtered_file ):
    """
    Filter an M5 alignment file to contain only the alignments with the fewest diffs
    """
    log.info('Filtering Blasr M5 results from "{0}"'.format( m5_file ))
    selected = {}
    diffs = {}
    count = 0
    for record in BlasrReader( m5_file ):
        count += 1
        diff_count = int(record.nmis) + int(record.nins) + int(record.ndel)
        if record.qname not in diffs:
            selected[record.qname] = record
            diffs[record.qname] = diff_count
        elif diffs[record.qname] > diff_count:
            results[record.qname] = record
            diffs[record.qname] = diff_count
    log.info('Selected %s records from %s alignments' % (count, len(selected)))
    with open( filtered_file, 'w' ) as output:
        for record in selected.itervalues():
            output.write('%s\n' % blasr_to_string( record ))
    log.info('Finished filtering Blasr results')

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
            filename, locus = line.strip().split()
            records += list( FastaReader( filename ) )
    return records
