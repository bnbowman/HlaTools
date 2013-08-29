import os

from pbhla.arguments import args
from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import fasta_size
from pbhla.utils import check_output_file, remove_file

NPROC = 6 

def align_best_reference( query, reference, output=None ):
    """
    Align the output of AA to the references and return
    """
    print args
    print args.nproc
    # Figure out the output and remove it if it exists
    if output is None:
        basename = '.'.join( query.split('.')[:-1] )
        output = '%s.m1' % basename
    remove_file( output )
    # Run Blasr
    ref_count = fasta_size( reference )
    blasr_args = {'nproc': 1,
                  'out': output,
                  'bestn': 1,
                  'nCandidates': ref_count,
                  'noSplitSubreads': True }
    run_blasr( query, reference, blasr_args )
    # Check the output file
    check_output_file( output )
    return output
