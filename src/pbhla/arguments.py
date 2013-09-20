import os, argparse

from . import __VERSION__
from pbhla.fasta.utils import is_fasta

# Default values for optiond
NUM_PROC = 8
MIN_SCORE = 0.75
MIN_LENGTH = 2000
MAX_COUNT = None

#EXCLUDED = ['DPB2', 'DRB3', 'DRB4', 'DRB5']
EXCLUDED = ['DPB2']
#SPLIT = ['DPB1', 'DRB1']
SPLIT = []

args = argparse.Namespace()

def parse_args():
    """
    Parse the options for running the HLA pipeline and
    """
    desc = "A pipeline for performing HLA haplotype sequencing."
    parser = argparse.ArgumentParser( description=desc )

    add = parser.add_argument
    add("input_file", 
        metavar="INPUT",
        help="A Fasta, BasH5 or FOFN of HLA data to haplotype.")
    add("reference_file", 
        metavar="REFERENCE", 
        help="FOFN of reference fasta files and their associated loci")
    add("genome", 
        metavar="GENOME", 
        help="Fasta file of the Human Genome")
    add("--raw_data",
        metavar="FILE",
        help='A BasH5 or FOFN of raw HLA data.  Only used if Input is Fasta.')
    add("--output", 
        metavar="DIR", 
        help="Destination folder for process results")
    add("--nproc", 
        metavar='INT', 
        type=int,
        default=NUM_PROC,
        help="Number of processors to use for parallelization ({0})".format(NUM_PROC))
    add("--min_read_score", 
        metavar='FLOAT', 
        type=float,
        default=MIN_SCORE,
        help="Only use reads with ReadScore greater than this ({0})".format(MIN_SCORE))
    add("--min_read_length", 
        metavar='INT', 
        type=int,
        default=MIN_LENGTH,
        help="Only use subreads longer than this ({0})".format(MIN_LENGTH))
    add("--max_count", 
        metavar='INT', 
        type=int,
        default=MAX_COUNT,
        help="Maximum number of subreads to use ({0})".format(MAX_COUNT))
    add("--smrt_path", 
        metavar="PATH", 
        help="Path to the setup script for the local SMRT Analysis installation")
    add("--resequence", 
        action="store_true", 
        help="Use quiver to resequence the top selected alleles")
    add("--amplicon_assembly",
        action="store_true",
        help="Use Amplicon Assembly to phase the Contigs")
    add("--exclude",
        nargs='*',
        metavar='LOCUS',
        default=EXCLUDED,
        help="List of loci to not include in the final summary or phase")
    add("--split",
        nargs='*',
        metavar='LOCUS',
        default=SPLIT,
        help="List of loci to split into two contigs")
    add("--msa", 
        metavar='FOFN',
        default=None,
        help="FOFN of prealigned MSAs and their associated loci")
    #add("--annotate", action='store_true', help="Avoid annotation.")

    class PrintVersionAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print "\tHLA Analysis Pipeline version: %s" % __VERSION__
            raise SystemExit

    add("--version",
        nargs=0,
        action=PrintVersionAction)

    parser.parse_args( namespace=args )

    args.input_file = os.path.abspath( args.input_file )

    if is_fasta( args.input_file ) and args.raw_data is None:
        msg = "Raw Data option must be specified in Input is Fasta"
        log.error( msg )
        raise ValueError( msg )
    elif args.raw_data:
        args.raw_data = os.path.abspath( args.raw_data )
    else:
        args.raw_data = args.input_file

    if args.output is None:
        args.output = args.input_file.split('.')[0]
