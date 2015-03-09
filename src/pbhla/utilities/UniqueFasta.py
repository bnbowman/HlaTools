#! /usr/bin/env python

import sys
from pbcore.io import FastaReader, FastaWriter, FastaRecord

seqs = set()
with FastaWriter( sys.stdout ) as handle:
    for rec in FastaReader( sys.argv[1] ):

        # Strip gaps from the sequence
        seq = rec.sequence.replace('-','').replace('.', '')

        # If we've seen the ungapped sequence before, skip
        if seq in seqs:
            continue
        # Otherwise add the sequence and write the record
        else:
            new_rec = FastaRecord(rec.name, seq)
            seqs.add( seq )
            handle.writeRecord( new_rec )
