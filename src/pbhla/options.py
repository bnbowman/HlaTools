import os, argparse

from . import __VERSION__
from pbhla.external.SmrtAnalysisTools import SmrtAnalysisRunner

# Default values for optiond
SMRT_ANALYSIS = "/mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh"
DILUTION = 1.0
MIN_SCORE = 0.8
MIN_LENGTH = 2500
NUM_PROC = 8

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
        help="A BasH5 or FOFN of BasH5s to haplotype.")
    add("reference_file", 
        metavar="REFERENCE", 
        help="FOFN of reference fasta files and their associated loci")
    add("genome", 
        metavar="GENOME", 
        help="Fasta file of the Human Genome")
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
    add("--dilution", 
        metavar='FLOAT', 
        type=float,
        default=DILUTION,
        help="Fraction of subreads to use ({0})".format(DILUTION))
    add("--smrt_path", 
        metavar="PATH", 
        default=SMRT_ANALYSIS,
        help="Path to the setup script for the local SMRT Analysis installation")
    add("--resequence", 
        action="store_true", 
        help="Use quiver to resequence the")
    #add("--MSA", help="FOFN of prealigned MSAs and their associated loci.")
    #add("--region_table", help="Region Table of White-Listed reads to use")
    #add("--phasr-args", nargs='*', default=[''], help="pass these args to phasr.")
    #add("--avoid_phasr", action='store_true',
    #    help="Avoid phasr if reference mapping has identified two likely phases.")
    #add("--annotate", action='store_true', help="Avoid annotation.")

    class PrintVersionAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print "    HLA Analysis Pipeline version: %s" % __VERSION__
            raise SystemExit

    add("--version",
        nargs=0,
        action=PrintVersionAction)

    # 
    parser.parse_args( namespace=args )

    args.input_file = os.path.abspath( args.input_file )
    if args.output is None:
        args.output = args.input_file.split('.')[0]

    # Check dilution factors
    if args.dilution <= 0.0 or args.dilution > 1.0:
        parser.error("Dilute factor must be between 0 and 1")

    # Initialize a tool for running SmrtAnalysis if needed
    if args.resequence:
        args.smrt_analysis = SmrtAnalysisRunner( self.smrt_path, self.log_files, args.nproc )

    # parse phasr args
    #self.phasr_argstring = ''
    #for argument in self.phasr_args:
    #    if ':' in argument:
    #        param, value = argument.split(":")
    #        self.phasr_argstring += '--%s %s ' % (param, value)
    #    elif 'output' in argument or 'cname' in argument:
    #        pass    
    #    else:
    #        self.phasr_argstring += '--%s ' % argument
