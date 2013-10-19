import os, argparse

from . import __VERSION__

NUM_PROC = 4

args = argparse.Namespace()

def parse_args():
    """
    Parse the options for running the HLA pipeline and
    """
    desc = "A pipeline for performing HLA haplotype sequencing."
    parser = argparse.ArgumentParser( description=desc )

    # Add the options and arguments
    add = parser.add_argument
    add("input_file", 
        metavar="INPUT",
        help="A Fasta, BasH5 or FOFN of HLA data to haplotype.")
    add("config_file",
        metavar="CONFIG",
        help="FOFN of reference fasta files and their associated loci")
    add("--output",
        metavar="DIR", 
        help="Destination folder for process results")
    add("--nproc", 
        metavar='INT', 
        type=int,
        default=NUM_PROC,
        help="Number of processors to use for parallelization ({0})".format(NUM_PROC))
    add("--version",
        nargs=0,
        action=PrintVersionAction)

    # Parse the options and update a few of the variables as needed
    parser.parse_args( namespace=args )

    args.input_file = os.path.abspath( args.input_file )

    if not args.output:
        args.output = args.input_file.split('.')[0]


class PrintVersionAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print "\tHLA Analysis Pipeline version: %s" % __VERSION__
        raise SystemExit