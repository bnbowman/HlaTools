#! /usr/bin/env python

import os, logging

from pbcore.io import FastaReader, FastqReader
from pbhla.filenames import get_file_type

from pbhla.io.BlasrIO import BlasrReader
from pbhla.cdna.extract_exons import extract_exons
from pbhla.cdna.exons_to_cDNA import exons_to_cDNA
from pbhla.sequences.utils import combine_sequences, write_sequences
from pbhla.utils import check_output_file, create_directory, remove_directory

log = logging.getLogger()

def extract_cDNA( input_file, exon_fofn, output=None,
                                         directory=None,
                                         reference_file=None, 
                                         alignment_file=None ):
    """
    Extract the cDNA sequences from a mixed Fasta
    """
    # Check the input files, and align the input file if needed
    if reference_file and alignment_file is None:
        alignment_file = align_best_reference( input_file, reference_file )
    elif reference_file is None and alignment_file is None:
        msg = "extract_alleles requires either an Alignment or a Reference!"
        log.error( msg )
        raise IOError( msg )

    # Set the output and directory if it hasn't been specified
    if directory is None:
        dirname = os.path.dirname( input_file )
        directory = os.path.join( dirname, 'cDNA' )
        remove_directory( directory )

    create_directory( directory )
    output = output or _get_output_file( input_file )

    # Prepare the Fasta by orienting and subsetting it
    records = _parse_input_records( input_file )
    fofn = _parse_exon_fofn( exon_fofn )
    loci = _parse_loci( alignment_file )
    log.info("Extracting cDNA sequences from all records")
    _extract_cDNA( records, loci, fofn, directory )
    log.info("Collecting all extracted cDNA records into %s" % output)
    _collect_cDNA( directory, output )

    # Clean up the directory and return the combined cDNA file
    remove_directory( directory )
    return output

def _get_output_file( input_file ):
    basename = '.'.join( input_file.split('.')[:-1] )
    file_type = get_file_type( input_file )
    return '%s.cDNA.%s' % (basename, file_type)

def _parse_exon_fofn( exon_fofn ):
    """
    Parse the Locus-specific FOFNs of Exons to a dictionary
    """
    fofn = {}
    with open( exon_fofn, 'r' ) as handle:
        for line in handle:
            locus, filename = line.strip().split()
            if locus in fofn:
                msg = 'Duplicate locus fofn "%s"!' % locus
                log.error( msg )
                raise ValueError( msg )
            else:
                fofn[locus] = filename
    return fofn

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

def _extract_cDNA( records, loci, fofn, directory ):
    """
    Extract and create a cDNA record for each Fasta Sequence
    """
    for record in records:
        # Create an output folder for each record to process
        name = record.name.split()[0]
        try:
            locus = loci[name]
        except:
            log.warn( 'No HLA locus associated with "%s" - skipping' % name )
            continue
        # Create a directory 
        record_directory = os.path.join( directory, name )
        create_directory( record_directory )
        # Find the appropriate Locus and FOFN
        if locus in fofn:
            exon_fofn = fofn[locus]
        else:
            log.warn( 'No exonic reference for %s' % locus )
        # Extract the exons and make the cDNA
        exon_file = extract_exons( record, exon_fofn, record_directory )
        if exon_file:
            cDNA_file = exons_to_cDNA( exon_file )

def _collect_cDNA( folder, output_file ):
    """
    Collect all of the cDNA sequences into one Fasta File
    """
    cDNA_files = []
    for entry in os.listdir( folder ):
        entry_path = os.path.join( folder, entry )
        if os.path.isdir( entry_path ):
            for filename in os.listdir( entry_path ):
                if filename.endswith('cDNA.fasta') or filename.endswith('cDNA.fastq'):
                    filepath = os.path.join( entry_path, filename )
                    cDNA_files.append( filepath )
    combine_sequences( cDNA_files, output_file )
    check_output_file( output_file )

if __name__ == '__main__':
    import sys
    logging.basicConfig( level=logging.INFO )

    input_file = sys.argv[1]
    exon_fofn = sys.argv[2]
    output = sys.argv[3]
    alignment_file = sys.argv[4]

    extract_cDNA( input_file, exon_fofn, output=output, 
                                         alignment_file=alignment_file )
