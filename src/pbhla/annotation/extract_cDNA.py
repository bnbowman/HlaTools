#! /usr/bin/env python

import sys, os, logging

from pbcore.io.FastaIO import FastaReader
from pbhla.io.BlasrIO import BlasrReader
from pbhla.fasta.update_orientation import update_orientation
from pbhla.fasta.utils import write_fasta, combine_fasta
from pbhla.utils import check_output_file, create_directory
from pbhla.annotation.extract_exons import extract_all_exons
from pbhla.annotation.exons_to_cDNA import exons_to_cDNA

log = logging.getLogger()

def extract_cDNA( fasta_file, blasr_file, exon_fofn, output="Typing" ):
    """
    Extract the cDNA sequences from a mixed Fasta
    """
    create_directory( output )
    input_fasta = _prepare_input_fasta( fasta_file, blasr_file, output )
    loci = _parse_loci( blasr_file )
    fofn = _parse_exon_fofn( exon_fofn )
    _extract_cDNA( input_fasta, loci, fofn, output )
    cDNA_file = _collect_cDNA( output )
    return cDNA_file

def _parse_exon_fofn( exon_fofn ):
    """
    Parse the Locus-specific FOFNs of Exons to a dictionary
    """
    fofn = {}
    with open( exon_fofn, 'r' ) as handle:
        for line in handle:
            filename, locus = line.strip().split()
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

def _prepare_input_fasta( fasta_file, blasr_file, output ):
    """
    Create an updated fasta with all sequences in the forward orientation
    """
    output_file = os.path.join( output, 'Oriented_Sequences.fasta' )
    update_orientation( fasta_file, blasr_file, output_file )
    check_output_file( output_file )
    return output_file

def _extract_cDNA( fasta_file, loci, fofn, output ):
    """
    Extract and create a cDNA record for each Fasta Sequence
    """
    for record in FastaReader( fasta_file ):
        # Create a directory for each sequence
        name = record.name.split()[0]
        record_dir = os.path.join( output, name )
        create_directory( record_dir )
        # Write out each sequence individually 
        query_file = "%s.fasta" % name
        query_path = os.path.join( record_dir, query_file )
        write_fasta( [record], query_path )
        # Find the appropriate Locus and FOFN
        locus = loci[name]
        if locus in fofn:
            exon_fofn = fofn[locus]
        else:
            continue
        # Extract the exons and make the cDNA
        exon_file = extract_all_exons( query_path, exon_fofn )
        cDNA_file = exons_to_cDNA( exon_file )

def _collect_cDNA( output ):
    """
    Collect all of the cDNA sequences into one Fasta File
    """
    cDNA_files = []
    for entry in os.listdir( output ):
        entry_path = os.path.join( output, entry )
        if os.path.isdir( entry_path ):
            for filename in os.listdir( entry_path ):
                if filename.endswith( '_cDNA.fasta' ):
                    filepath = os.path.join( entry_path, filename )
                    cDNA_files.append( filepath )
    output_file = os.path.join( output, "cDNA_Sequences.fasta" )
    combine_fasta( cDNA_files, output_file )
    check_output_file( output_file )
    return output_file

if __name__ == '__main__':
    logging.basicConfig( level=logging.INFO )

    fasta_file = sys.argv[1]
    blasr_file = sys.argv[2]
    exon_fofn = sys.argv[3]
    output = sys.argv[4]

    extract_cDNA( fasta_file, blasr_file, exon_fofn, output ) 
