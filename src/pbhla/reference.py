import os, re, logging

from pbcore.io.FastaIO import FastaReader
from pbhla.fasta.utils import write_fasta

log = logging.getLogger()

def parse_reference_fofn( fofn_file ):
    basename = os.path.basename( fofn_file )
    log.info('Parsing the reference FOFN: "%s"' % basename)
    sequences = _parse_reference_sequences( fofn_file )
    metadata = _parse_reference_metadata( fofn_file )
    loci = _parse_reference_loci( fofn_file )
    log.info('Finished parsing the reference FOFN')
    return sequences, metadata, loci

def _parse_reference_sequences( fofn_file ):
    log.info('Parsing reference sequence data...')
    records = []
    with open( fofn_file, 'r') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            filename, locus = line.strip().split()
            records += list( FastaReader( filename ) )
    log.info("Found %s reference sequence records" % len(records))
    return records

def _parse_reference_metadata( fofn_file ):
    log.info('Parsing reference metadata...')
    metadata = {}
    with open( fofn_file, 'r') as handle:
        for line in handle:
            if not line.startswith('#'):
                continue
            key, value = line[1:].strip().split('=')
            metadata[key] = value
    log.info('Found %s pieces of metadata' % len(metadata))
    return metadata

def _parse_reference_loci( fofn_file ):
    log.info('Parsing reference loci...')
    results = {}
    loci = set()
    with open(fofn_file, 'r') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fasta_path, locus = line.strip().split()
            loci.add( locus )
            fasta_file = os.path.basename( fasta_path )
            log.info('Reading "{0}" sequences from "{1}"'.format(locus, fasta_file))
            for record in FastaReader( fasta_path ):
                name = record.name.split()[0]
                if not re.search('__', name):
                    name = name.split('_')[0]
                results[name] = locus
    log.info('Finished parsing data for {0} sequences from {1} loci'.format( len(results), 
                                                                             len(loci) ))
    return results
