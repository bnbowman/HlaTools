#! /usr/bin/env python

from pbcore.io import FastaReader, FastaWriter, FastaRecord

def split_imgt_drb( input_file ):
    drb3 = FastaWriter("DRB3_nuc.fasta")
    drb4 = FastaWriter("DRB4_nuc.fasta")
    drb5 = FastaWriter("DRB5_nuc.fasta")

    for record in FastaReader( input_file ):
        # Check that this is an IMGT-formatted FASTA record
        assert record.header.startswith('HLA:')

        # Split the locus name out of the header
        locus = record.header.strip().split('_')[1]

        # Write the record to the appropriate file
        if locus.startswith('DRB3'):
            drb3.writeRecord( record )
        elif locus.startswith('DRB4'):
            drb4.writeRecord( record )
        elif locus.startswith('DRB5'):
            drb5.writeRecord( record )
        else:
            raise Exception("Locus {0} is not DRB3, 4, or 5!".format(locus))

if __name__ == '__main__':
    import sys

    input_file = sys.argv[1]

    split_imgt_drb( input_file )
