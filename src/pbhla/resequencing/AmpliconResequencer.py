#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
from pbcore.io import FastaReader
from pbhla.resequencing.Resequencer import Resequencer
from pbhla.fasta.utils import is_fasta
from pbhla.utils import get_file_root

def AmpliconResequencer( amplicon_dir, barcoded_dir, output="AmpliconResequencer", setup=None, nproc=1 ):
    alleles = get_consensus_alleles( amplicon_dir )
    barcodes = sort_alleles_by_barcode( alleles )
    barcoded_read_files = get_barcoded_read_files( barcoded_dir, barcodes.keys() )
    for b, a in barcoded_read_files.iteritems():
        print b, a

    resequencer = Resequencer( setup=setup, nproc=nproc )



def get_consensus_alleles( amplicon_folder ):
    # Parse the 'Good' or Non-Chimeric alleles
    good_output = os.path.join( amplicon_folder, 'amplicon_analysis.fasta' )
    good_alleles = [r for r in FastaReader( good_output )]
    good_names = set([r.name.split()[0] for r in good_alleles])
    assert len(good_names) == len(good_alleles)

    # Parse the 'Bad' or Chimeric/Noise alleles
    bad_output = os.path.join( amplicon_folder, 'amplicon_analysis_chimeras_noise.fasta')
    bad_alleles = [r for r in FastaReader( bad_output )]
    bad_names = set([r.name.split()[0] for r in bad_alleles])
    assert len(bad_names) == len(bad_alleles)

    # Make sure the two sets don't overlap and return the results
    assert len(good_names & bad_names) == 0
    combined_alleles = good_alleles + bad_alleles
    return combined_alleles

def sort_alleles_by_barcode( alleles ):
    def get_barcode( record ):
        name = record.name.split()[0]
        return name.split('_')[0][7:]

    barcodes = {}
    for record in alleles:
        barcode = get_barcode( record )
        try:
            barcodes[barcode].append( record )
        except:
            barcodes[barcode] = [record]
    return barcodes

def get_barcoded_read_files( barcoded_dir, barcodes ):
    barcoded_read_files = {}
    for filename in os.listdir( barcoded_dir ):
        if is_fasta( filename ):
            rootname = get_file_root( filename )
            if rootname in barcodes:
                filepath = os.path.join( barcoded_dir, filename )
                barcoded_read_files[rootname] = filepath
    return barcoded_read_files


if __name__ == '__main__':
    import sys

    amplicon_dir = sys.argv[1]
    barcoded_dir = sys.argv[2]
    output = sys.argv[3]
    setup = sys.argv[4]
    nproc = int(sys.argv[5])

    AmpliconResequencer( amplicon_dir, barcoded_dir, output=output,
                                                     setup=setup,
                                                     nproc=nproc )