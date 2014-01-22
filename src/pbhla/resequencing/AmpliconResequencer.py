#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
from pbcore.io import FastaReader
from pbhla.resequencing.Resequencer import Resequencer

def AmpliconResequencer( amplicon_folder, output="AmpliconResequencer", setup=None, nproc=1 ):
    alleles = get_consensus_alleles( amplicon_folder )
    print len(alleles)

    resequencer = Resequencer( setup=setup, nproc=nproc )

def get_consensus_alleles( amplicon_folder ):
    # Parse the 'Good' or Non-Chimeric alleles
    good_output = os.path.join( amplicon_folder, 'amplicon_analysis.fasta' )
    good_alleles = [r for r in FastaReader( good_output )]
    good_names = set([r.name.split() for r in good_alleles])
    assert len(good_names) == len(good_alleles)

    # Parse the 'Bad' or Chimeric/Noise alleles
    bad_output = os.path.join( amplicon_folder, 'amplicon_analysis_chimeras_noise.fasta')
    bad_alleles = [r for r in FastaReader( bad_output )]
    bad_names = set([r.name.split() for r in bad_alleles])
    assert len(bad_names) == len(bad_alleles)

    # Make sure the two sets don't overlap and return the results
    assert len(good_names & bad_names) == 0
    combined_alleles = good_alleles + bad_alleles
    return combined_alleles

if __name__ == '__main__':
    import sys

    amplicon_folder = sys.argv[1]
    output_folder = sys.argv[2]
    setup = sys.argv[3]
    nproc = int(sys.argv[4])

    AmpliconResequencer( amplicon_folder, output=output,
                                          setup=setup,
                                          nproc=nproc )