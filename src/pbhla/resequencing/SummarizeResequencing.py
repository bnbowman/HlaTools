#! /usr/bin/env python

import os, re, sys, logging

from shutil import copy

from pbhla.utils import create_directory, write_list_file
from pbhla.fasta.utils import fasta_size, fasta_length, copy_fasta

log = logging.getLogger()

def combine_resequencing_output(input_dir, output_dir):
    assert os.path.isdir( input_dir )
    create_directory( output_dir )
    log.info('Combining resequencing output from "%s" in "%s"' % (input_dir, output_dir))

    resequencing_results = find_resequencing_results( input_dir )
    fasta_files, fastq_files = output_results( resequencing_results, output_dir )

    fasta_fofn = os.path.join( output_dir, 'Resequenced_Clusense_Fasta.fofn' )
    write_list_file( fasta_files, fasta_fofn )

    fastq_fofn = os.path.join( output_dir, 'Resequenced_Clusense_Fastq.fofn' )
    write_list_file( fastq_files, fastq_fofn )
    return fasta_fofn, fastq_fofn

def find_resequencing_results( input_dir ):
    log.info('Identifying individual Resequencing output folders in "{0}"'.format(input_dir))
    resequencing_results = []
    for entry in os.listdir( input_dir ):
        entry_path = os.path.join( input_dir, entry )
        fasta, fastq = get_consensus_seqs( entry_path )
        if fasta and fastq:
            resequencing_results.append( (entry, fasta, fastq) )
    log.info('Identified %s successfully resequenced contigs' % len(resequencing_results))
    return resequencing_results

def get_consensus_seqs( dir_name ):
    if not os.path.isdir( dir_name ):
        return ( False, False )
    fasta_path = os.path.join( dir_name, 'consensus.fasta' )
    fastq_path = os.path.join( dir_name, 'consensus.fastq' )
    if os.path.isfile( fasta_path ) and os.path.isfile( fastq_path ):
        return ( fasta_path, fastq_path )
    return ( False, False )

def output_results( results, output_dir ):
    log.info('Outputting resequencing results to "{0}"'.format(output_dir))
    fasta_files = []
    fastq_files = []
    for contig, fasta, fastq in results:
        # Copy the Fasta consensus
        output_fasta = contig + '.fasta'
        fasta_path = os.path.join( output_dir, output_fasta )
        copy(fasta, fasta_path)
        fasta_files.append( fasta_path )
        # Copy the Fastq consensus
        output_fastq = contig + '.fastq'
        fastq_path = os.path.join( output_dir, output_fastq )
        copy(fastq, fastq_path)
        fastq_files.append( fastq_path )
    log.info('Finished Outputting resequencing results')
    return (fasta_files, fastq_files)

if __name__ == "__main__":
    import argparse
    desc = "Combine a series of Resequencing directories"
    parser = argparse.ArgumentParser(description=desc)

    add = parser.add_argument
    add("input_dir", metavar="INPUT_DIR", 
        help="Fasta-format file of the reference sequence to use")
    add("output_dir", metavar="OUTPUT_DIR", 
        help="Name of the directory to output results to")
    args = parser.parse_args()

    combine_clusense_output( args.input_dir, args.output_dir )
