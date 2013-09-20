import os

from pbhla.arguments import args
from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import fasta_size
from pbhla.utils import check_output_file, remove_file

if hasattr( args, 'nproc' ):
    NPROC = args.nproc
else:
    NPROC = 6

def align_best_reference( query, reference, output=None ):
    """
    Align the output of AA to the references and return
    """
    # Figure out the output and remove it if it exists
    if output is None:
        basename = '.'.join( query.split('.')[:-1] )
        output = '%s.m1' % basename
    remove_file( output )
    # Run Blasr
    ref_count = fasta_size( reference )
    blasr_args = {'nproc': NPROC,
                  'out': output,
                  'bestn': 1,
                  'nCandidates': ref_count,
                  'noSplitSubreads': True }
    run_blasr( query, reference, blasr_args )
    # Check the output file
    check_output_file( output )
    return output

def full_align_best_reference( query, reference, output=None ):
    """
    Align the output of AA to the references and return
    """
    # Figure out the output and remove it if it exists
    if output is None:
        basename = '.'.join( query.split('.')[:-1] )
        output = '%s.m5' % basename
    remove_file( output )
    # Run Blasr
    ref_count = fasta_size( reference )
    blasr_args = {'nproc': NPROC,
                  'out': output,
                  'm': 5,
                  'bestn': 1,
                  'nCandidates': ref_count,
                  'noSplitSubreads': True }
    run_blasr( query, reference, blasr_args )
    # Check the output file
    check_output_file( output )
    return output
