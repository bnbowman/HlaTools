#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

#Test
from pbhla.resequencing.AmpAnalysisResequencer import AmpliconAnalysisResequencer

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    add = parser.add_argument
    add('data_file', metavar='DATA',
        help="BasH5 or FOFN of sequence data")
    add('barcode_file', metavar='BARCODE',
        help="BasH5 or FOFN of sequence data")
    add('amplicon_analysis', metavar='INPUT',
        help="BasH5 or FOFN of sequence data")
    add('data_file', metavar='DATA',
        help="BasH5 or FOFN of sequence data")
    add('--output', default='resequencing', metavar='DIR',
        help="Specify a directory for intermediate files")
    add('--setup', metavar='SETUP_FILE',
        help='Path to the SMRT Analysis setup script')
    add('--nproc', type=int, default=1, metavar='INT',
        help="Number of processors to use")
    args = parser.parse_args()

    # Run the specified resequencing process
    resequencer = AmpliconAnalysisResequencer( setup=args.setup,
                                               nproc=args.nproc )
    resequencer( args.data_file,
                 args.barcode_file,
                 args.amplicon_analysis,
                 output=args.output,
                 barcodes=args.barcode_list )