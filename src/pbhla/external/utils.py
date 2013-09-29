import os

from pbhla.arguments import args, NUM_PROC
from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import fasta_size
from pbhla.utils import check_output_file, remove_file, valid_file

if hasattr( args, 'nproc' ):
    nproc = args.nproc
else:
    nproc = NUM_PROC

def get_alignment_file( query, reference_file, alignment_file ):
    if alignment_file and valid_file( alignment_file ):
        return alignment_file
    elif reference_file:
        return align_best_reference( query, reference_file )
    msg = "No valid Alignment or a Reference files found"
    log.error( msg )
    raise IOError( msg )

def align_best_reference( query, reference, output=None ):
    """
    Align the output of AA to the references and return
    """
    output = output or _get_output_file( query, 'm1' )
    remove_file( output )
    # Run Blasr
    ref_count = fasta_size( reference )
    blasr_args = {'nproc': nproc,
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
    output = output or _get_output_file( query, 'm5' )
    remove_file( output )
    # Run Blasr
    ref_count = fasta_size( reference )
    blasr_args = {'nproc': nproc,
                  'out': output,
                  'm': 5,
                  'bestn': 1,
                  'nCandidates': ref_count,
                  'noSplitSubreads': True }
    run_blasr( query, reference, blasr_args )
    # Check the output file
    check_output_file( output )
    return output

def _get_output_file( query, suffix ):
    """
    Get an output filename of the appropriate type
    """
    basename = '.'.join( query.split('.')[:-1] )
    output = '%s.%s' % (basename, suffix)
    return output
