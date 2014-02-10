#! /usr/bin/env python

from pbhla.typing.sequences import type_sequences

if __name__ == '__main__':
    import sys

    input = sys.argv[1]
    grouping = sys.argv[2] if len(sys.argv) > 2 else None
    exon_fofn = sys.argv[3] if len(sys.argv) > 3 else None
    genomic_reference = sys.argv[4] if len(sys.argv) > 4 else None
    cDNA_reference = sys.argv[5] if len(sys.argv) > 5 else None

    type_sequences( input, grouping, exon_fofn, genomic_reference, cDNA_reference )
