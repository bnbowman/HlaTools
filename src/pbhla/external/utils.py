from pbhla.arguments import args
from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import fasta_size
from pbhla.utils import check_output_file, remove_file

log = logging.getLogger(__name__)

if hasattr(args, 'nproc'):
    NPROC = args.nproc
else:
    NPROC = 6


def _get_output_file(query, output, suffix):
    """
    Get the output file for a call to Blasr
    """
    if output is None:
        basename = '.'.join(query.split('.')[:-1])
        output = '%s.%s' % (basename, suffix)
    remove_file(output)
    return output


def align_best_reference(query, reference, output=None):
    """
    Align the output of AA to the references and return
    """
    output = _get_output_file(query, output, '.m1')
    # Run Blasr
    ref_count = fasta_size(reference)
    blasr_args = {'nproc': NPROC,
                  'out': output,
                  'bestn': 1,
                  'nCandidates': ref_count,
                  'noSplitSubreads': True}
    run_blasr(query, reference, blasr_args)
    # Check the output file
    check_output_file(output)
    return output


def full_align_best_reference(query, reference, output=None):
    """
    Align the output of AA to the references and return
    """
    # Figure out the output and remove it if it exists
    output = _get_output_file(query, output, '.m5')
    # Run Blasr
    ref_count = fasta_size(reference)
    blasr_args = {'nproc': NPROC,
                  'out': output,
                  'm': 5,
                  'bestn': 1,
                  'nCandidates': ref_count,
                  'noSplitSubreads': True}
    run_blasr(query, reference, blasr_args)
    # Check the output file
    check_output_file(output)
    return output
