#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'


from pbhla.resequencing.AmpAnalysisResequencer import AmpliconAnalysisResequencer

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    add = parser.add_argument
    add('amplicon_analysis', metavar='INPUT',
        help="Fasta/Fastq/Folder of Amplicon Analysis output")
    add('data_file', metavar='DATA',
        help="BasH5 or FOFN of sequence data")
    add('barcode_file', metavar='BARCODE',
        help="BcH5 or FOFN of barcode data")
    add('-s', '--setup', metavar='SETUP_FILE',
        help='Path to the SMRT Analysis setup script')
    add('-l', '--min_length', metavar='INT',
        type=int, default=500,
        help='Minimum subread length to use for resequencing')
    add('-r', '--min_snr', metavar='FLOAT',
        type=float, default=4.0,
        help='Minimum Signal-To-Noise Ratio to use for resequencing')
    add('-b', '--barcode_list', metavar='BARCODE_LIST',
        help="Comma-separated list of barcodes to resequence")
    add('-o', '--output', default='resequencing', metavar='DIR',
        help="Specify a directory for intermediate files")
    add('--nproc', type=int, default=1, metavar='INT',
        help="Number of processors to use")
    add('--debug', action='store_true',
        help="Flag to enable Debug mode")
    args = parser.parse_args()

    # Run the specified resequencing process
    resequencer = AmpliconAnalysisResequencer( args.output,
                                               setup=args.setup,
                                               nproc=args.nproc,
                                               debug=args.debug )
    resequencer( args.amplicon_analysis,
                 args.data_file,
                 args.barcode_file,
                 barcode_string=args.barcode_list,
                 min_length=args.min_length,
                 min_snr=args.min_snr )