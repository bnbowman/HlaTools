#! /usr/bin/env python

import os, re, sys, logging

from shutil import copy

from pbhla.utils import create_directory, write_fofn
from pbhla.fasta.utils import fasta_size, fasta_length, copy_fasta

log = logging.getLogger()

CNS_FOFN = 'Clusense_Consensus_Files.txt'
READ_FOFN = 'Clusense_Read_Files.txt'

def combine_clusense_output(input_dir, output_dir):
    assert os.path.isdir( input_dir )
    create_directory( output_dir )
    log.info('Combining clusense output from "{0}" in "{1}"'.format(input_dir, output_dir))

    clusense_dirs = find_clusense_dirs( input_dir )
    clusense_clusters = find_clusense_clusters( clusense_dirs )
    cns_files, read_files = output_clusters( clusense_clusters, output_dir )

    cns_output = os.path.join( output_dir, CNS_FOFN )
    write_fofn( cns_files, cns_output )

    read_output = os.path.join( output_dir, READ_FOFN )
    write_fofn( read_files, read_output )
    return cns_output, read_output

def find_clusense_dirs( input_dir ):
    log.info('Identifying individual Clusense output folders in "{0}"'.format(input_dir))
    clusense_dirs = []
    for entry in os.listdir( input_dir ):
        entry_path = os.path.join( input_dir, entry )
        if is_clusense_dir( entry_path ):
            clusense_dirs.append( entry_path )
    if clusense_dirs:
        clusense_dirs = list(set( clusense_dirs ))
        log.info('Identified {0} Clusense output folders'.format(len(clusense_dirs)))
        return clusense_dirs
    elif is_clusense_dir( input_dir ):
        log.info('Identified 1 Clusense output folder')
        return [ input_dir ]
    log.info('Identified 0 Clusense output folders!')
    return []

def find_clusense_clusters( clusense_dirs ):
    log.info('Identifying individual Clusense clusters from identified output folders')
    clusense_clusters = []
    for entry in clusense_dirs:
        cluster = 1
        while True:
            result = cluster_exists( entry, cluster )
            if result:
                clusense_clusters.append( result )
            else:
                break
            cluster += 1
    log.info('Identified {0} individual Clusense clusters'.format(len(clusense_clusters)))
    return clusense_clusters

def output_clusters( clusters, output_dir ):
    log.info('Outputting high-quality clusters to "{0}"'.format(output_dir))
    cns_files = []
    read_files = []
    for consensus, reads in clusters:
        contig_dir = os.path.dirname( consensus )
        contig_name = os.path.split( contig_dir )[1]
        cluster = consensus.split('_')[-2]
        # Copy the consensus file
        cns_name = "{0}_{1}_cns".format(contig_name, cluster)
        cns_file = cns_name + '.fasta'
        cns_path = os.path.join( output_dir, cns_file )
        copy_fasta(consensus, cns_path, cns_name)
        cns_files.append( cns_path )
        # Copy the reads file
        read_file = "{0}_{1}.fasta".format(contig_name, cluster)
        read_path = os.path.join( output_dir, read_file )
        copy(reads, read_path)
        read_files.append( read_path )
    log.info('Finished Outputting high-quality clusters')
    return (cns_files, read_files)

def output_list( file_list, output_dir, filename ):
    log.info('Outputting file-list to "{0}"'.format(filename))
    output_path = os.path.join( output_dir, filename )
    with open(output_path, 'w') as handle:
        for file_path in file_list:
            print >> handle, file_path

def is_clusense_dir( dir_name ):
    if os.path.isdir( dir_name ):
        return cluster_exists( dir_name, 1 )
    return False

def cluster_exists( dir_name, cluster ):
    cluster = str(cluster).zfill(2)
    cns_file = 'group_' + cluster + '_cns.fa'
    cns_path = os.path.join( dir_name, cns_file )
    read_file = 'group_' + cluster + '.fa'
    read_path = os.path.join( dir_name, read_file )
    if os.path.isfile( cns_path ) and os.path.isfile( read_path ):
        return ( cns_path, read_path )
    return False

if __name__ == "__main__":
    import argparse
    desc = "Combine the "
    parser = argparse.ArgumentParser(description=desc)

    add = parser.add_argument
    add("input_dir", metavar="INPUT_DIR", 
        help="Fasta-format file of the reference sequence to use")
    add("output_dir", metavar="OUTPUT_DIR", 
        help="Name of the directory to output results to")
    args = parser.parse_args()

    combine_clusense_output( args.input_dir, args.output_dir )
