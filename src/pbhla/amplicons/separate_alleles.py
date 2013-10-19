#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

import logging

from collections import Counter

from pbcore.io.FastaIO import FastaReader, FastaWriter

from pbhla.external.utils import align_best_reference
from pbhla.io.BlasrIO import BlasrReader
from pbhla.utils import valid_file
from pbhla.typing.HlaTypes import HlaType

log = logging.getLogger(__name__)

def separate_alleles( subread_file, reference_file ):
    """
    Separate subreads from highly divergent alleles
    """
    ordered_refs = order_references( subread_file, reference_file )
    refs = select_references( reference_file, ordered_refs )
    reference_subset = subset_references( reference_file, refs )
    write_references( reference_file, refs )
    subread_refs = sort_subreads( subread_file, reference_subset )
    separate_subreads( subread_file, refs, subread_refs )

def order_references( subread_file, reference_file ):
    """
    Select the two best reference sequences from a list
    """
    log.info("Selecting the best references sequences to use")
    temp = 'temp.m1'
    if not valid_file( temp ):
        align_best_reference( subread_file, reference_file, temp )
    c = Counter([hit.tname for hit in BlasrReader(temp)])
    return [k for k, v in c.most_common()]

def select_references( reference_file, refs ):
    for i, ref in enumerate(refs):
        for record in FastaReader( reference_file ):
            if record.name.startswith( ref ):
                hla_type = HlaType.from_string(record.name)
                if i == 0:
                    first = ref
                    first_type = hla_type
                elif first_type.field1 != hla_type.field1:
                    return (first, ref)


def subset_references( reference_file, reference_names ):
    output = 'references.fasta'
    with FastaWriter(output) as writer:
        for record in FastaReader(reference_file):
            name = record.name.split()[0]
            if name in reference_names:
                writer.writeRecord( record )
    return output

def write_references( reference_file, references ):
    for i, ref in enumerate(references):
        for record in FastaReader(reference_file):
            name = record.name.split()[0]
            if name == ref:
                filename = 'reference_%s.fasta' % (i+1)
                with FastaWriter( filename ) as writer:
                    writer.writeRecord( record )

def sort_subreads( subread_file, reference_file ):
    """
    Aligning
    """
    log.info("Aligning subreads to the two best references")
    temp = 'temp2.m1'
    if valid_file( temp ):
        return {hit.qname: hit.tname for hit in BlasrReader(temp)}
    align_best_reference( subread_file, reference_file, temp )
    return {hit.qname: hit.tname for hit in BlasrReader(temp)}

def separate_subreads( subread_file, refs, subread_refs ):
    """
    Pass
    """
    writers = open_allele_writers( refs )
    for record in FastaReader( subread_file ):
        name = record.name.split()[0]
        try:
            ref = subread_refs[name]
        except KeyError:
            log.warn('"No reference found for "%s"' % name)
            continue
        writers[ref].writeRecord( record )
    close_allele_writers( writers )

def open_allele_writers( references ):
    writers = {}
    for i, ref in enumerate( references ):
        filename = 'allele_%s.fasta' % (i+1)
        writers[ref] = FastaWriter(filename)
    return writers

def close_allele_writers( writers ):
    for reference, writer in writers.iteritems():
        writer.close()

if __name__ == '__main__':
    import sys

    subread_file = sys.argv[1]
    reference_file = sys.argv[2]

    logging.basicConfig( stream=sys.stdout )
    separate_alleles( subread_file, reference_file )