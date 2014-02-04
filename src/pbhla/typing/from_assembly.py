#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

from pbcore.io import FastaRecord

from pbhla.fasta.utils import read_fasta_dict, reverse_complement, write_fasta
from pbhla.external.commandline_tools import run_blasr
from pbhla.utils import check_output_file
from pbhla.io.BlasrIO import BlasrReader

PLOIDY = 'haploid'

def from_assembly( contig_file, reference_fofn ):
    contigs = read_fasta_dict( contig_file )
    references = read_reference_fofn( reference_fofn )
    all_genes = []
    for locus, reference in references.iteritems():
        alignment = align_reference_to_contigs( locus, reference, contig_file )
        hits = read_blasr_hits( alignment )
        hits = pick_blasr_hits( hits )
        genes = extract_genes( contigs, hits )
        all_genes += genes
    write_fasta( all_genes, 'output.fasta' )

def read_reference_fofn( reference_fofn ):
    references = {}
    for line in open( reference_fofn ):
        if line.startswith('#'):
            continue
        reference, locus = line.strip().split()
        references[locus] = reference
    return references

def align_reference_to_contigs( locus, reference, contig_file ):
    output_file = 'HLA-%s.m1' % locus
    blasr_args = {'nproc': 8,
                  'out': output_file,
                  'bestn': 1,
                  'noSplitSubreads': True}
    run_blasr( reference, contig_file, blasr_args )
    check_output_file( output_file )
    return output_file

def read_blasr_hits( alignment ):
    hits = list(BlasrReader( alignment ))
    return sorted(hits, key=lambda x: float(x.pctsimilarity), reverse=True)

def pick_blasr_hits( hits ):
    hq_hits = [h for h in hits if float(h.pctsimilarity) > 98.0]
    if len(hq_hits) < 2:
        return [ hits[0] ]
    print [(h.tstart, h.tend) for h in hq_hits]
    long_hits = sorted(hq_hits, key=lambda x: hit_length(x), reverse=True)
    return [ long_hits[0] ]

def hit_length( hit ):
    return abs(int(hit.tend) - int(hit.tstart))

def extract_genes( contigs, hits ):
    genes = []
    for hit in hits:
        start = int(hit.tstart)
        end   = int(hit.tend)
        name  = '%s_%s_%s_NumReads10' % (hit.tname, hit.tstart, hit.tend)
        if hit.qstrand == hit.tstrand:
            contig = contigs[hit.tname]
        else:
            contig = reverse_complement( contigs[hit.tname] )
        sequence = contig.sequence[start:end+1]
        gene = FastaRecord(name, sequence)
        genes.append( gene )
    return genes


if __name__ == '__main__':
    import sys

    fasta_file = sys.argv[1]
    reference_fofn = sys.argv[2]
    ploidy = sys.argv[3] if len(sys.argv) > 3 else PLOIDY

    from_assembly( fasta_file, reference_fofn )