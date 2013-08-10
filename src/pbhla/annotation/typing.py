#! /usr/bin/env python

import os, re, sys, argparse, logging, subprocess

from collections import namedtuple

from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord
from pbcore.io.GffIO import GffReader, Gff3Record, GffWriter

from pbhla.io.GffIO import ( create_annotation,
                             create_annotation2,
                             create_var_annotation )
from pbhla.external.commandline_tools import run_muscle
from pbhla.fasta.utils import fasta_size
from pbhla.utils import ( create_directory, 
                          read_dict_file,
                          memoize )
from pbhla.fasta.utils import write_fasta
from pbhla.align.MultiSequenceAligner import MSA_aligner

log = logging.getLogger()

info = namedtuple('info', 'canonical_pos, feature, codon')

def type_hla( MSA_fofn, fasta_file, locus_dict, output_dir ):
    log.info("Initializing Annotation Process")

    create_directory( output_dir )

    if isinstance(locus_dict, str):
        locus_dict = read_dict_file( locus_dict )

    log.info("Found %s sequences to annotate" % len(locus_dict))

    MSA_fn_dict={} 
    MSA_cDNA_fn_dict={}
    MSA_info_fn_dict={} 
    MSA_cDNA_info_fn_dict={}

    temp_fasta = os.path.join( output_dir, 'tmp.fasta' )
    with open( MSA_fofn, "r") as f:
        for line in f:
            genomic, nucleotide, locus = line.split()
            MSA_fn_dict[locus] = genomic + ".afa"
            MSA_info_fn_dict[locus] = genomic + ".info"
            MSA_cDNA_info_fn_dict[locus] = nucleotide + ".info"
            MSA_cDNA_fn_dict[locus] = nucleotide + ".afa"

    gDNA_var_path = os.path.join(output_dir, 'best_gDNA_ref_match_variants.gff')
    gDNA_var_writer = GffWriter( gDNA_var_path )
    gDNA_var_writer.writeHeader('##pacbio-variant-version 1.4')

    cDNA_var_path = os.path.join(output_dir, 'best_cDNA_ref_match_variants.gff')
    cDNA_var_writer = GffWriter( cDNA_var_path )
    cDNA_var_writer.writeHeader('##pacbio-variant-version 1.4')

    gDNA_annot_path = os.path.join(output_dir, "gDNA.gff")
    gDNA_annot_writer = GffWriter( gDNA_annot_path )

    cDNA_annot_path = os.path.join(output_dir, "cDNA.gff")
    cDNA_annot_writer = GffWriter( cDNA_annot_path )

    feature_path = os.path.join(output_dir, "features.gff")
    feature_writer = GffWriter( feature_path )

    gDNA_allele_path = os.path.join(output_dir, 'gDNA_allele_calls.txt')
    gDNA_allele_handle = open(gDNA_allele_path, 'w')

    cDNA_allele_path = os.path.join(output_dir, 'cDNA_allele_calls.txt')
    cDNA_allele_handle = open(cDNA_allele_path, 'w')

    cDNA_path = os.path.join(output_dir, 'cDNA_sequences.fasta')
    cDNA_writer = FastaWriter( cDNA_path )

    for r in FastaReader( fasta_file ):
        print r
        write_fasta([r], temp_fasta)

        name = r.name.strip()
        print r.name, name
        if name.endswith('quiver'):
            name = name.split("|")[0]
        if name.endswith('_cns'):
            name = name[:-4]

        print locus_dict
        try:
            locus = locus_dict[name]
        except:
            continue

        if locus not in ['A', 'B', 'C', 'DPA', 'DQA', 'DQB']:
            continue
        log.info('Processing "%s" from locus "%s"' % (name, locus))

        ### read in profile features
        MSA_info_fn = MSA_info_fn_dict[locus]
        MSA_info = parse_genomic_info( MSA_info_fn )

        ### read profile fn
        MSA_fn = MSA_fn_dict[locus]

        ### align sequence to profile for this locus
        output = os.path.join(output_dir, name + ".afa")
        reference = get_reference_record( MSA_fn )
        add_to_alignment(temp_fasta, MSA_fn, output)
        aligned_reference = get_aligned_record( output, reference.name )
        aligned_query = get_aligned_record( output, r.name )

        ### create annotation gff3 for this sequence
        ### comparing the same ref between the new and old profile allows us to keep track of the "canonical" coordinates
        #gff3_annotation = create_annotation(comparison_seq, new_comparison_seq, consensus_seq, coverage_map, MSA_info, r.name, quality_map)     
        gff3_annotation = create_annotation2( reference.sequence,
                                              aligned_reference.sequence,
                                              aligned_query.sequence,
                                              MSA_info,
                                              r.name )
        last_feature = None
        feature_regions = {}
        for record in gff3_annotation:
            gDNA_annot_writer.writeRecord(record)
            try:
                feature_regions[record.feature].append( record.sequence_position )
            except KeyError:
                feature_regions[record.feature] = [int(record.sequence_position)]

        for item in feature_regions.iteritems():
            feature = item[0]
            if 'exon' not in feature: 
                continue
            bases = sorted(item[1])     
            start = bases[0]
            end = bases[-1]
            record = Gff3Record(r.name, start, end, 'region') 
            record.feature = feature
            feature_writer.writeRecord(record)
            
        ### find best match among all profiles for our consensus sequence       
        ### write it out to a file
        best_gDNA, match, total, score, gDNA_var_map = MSA_aligner( aligned_query.sequence, 
                                                                    output, 
                                                                    r.name )           
        gDNA_allele_handle.write("%s %s %s %s %s %s\n" % ( locus,
                                                           r.name, 
                                                           best_gDNA, 
                                                           total-match, 
                                                           total, 
                                                           100*float(score) ))
        
        ### create artificial cDNA sequence from the gDNA sequence
        ### create a dict between gDNA and cDNA position so that we can read in quality and coverage info to the cDNA gff
        ### and so we can keep track of where our artificial cDNA sequence fits in to the real sequence
        cDNA_consensus_sequence = ''
        cDNA_base_counter = 1
        gDNA_pos_to_cDNA_pos = {}

        for i in xrange(len(gff3_annotation)):
            if gff3_annotation[i].type == 'base':       
                if gff3_annotation[i].isexon:
                    cDNA_consensus_sequence += gff3_annotation[i].basecall
                    gDNA_pos_to_cDNA_pos[gff3_annotation[i].sequence_position] = cDNA_base_counter
                    cDNA_base_counter += 1

        cDNA_pos_to_gDNA_pos={}
        for item in gDNA_pos_to_cDNA_pos.iteritems():
            cDNA_pos_to_gDNA_pos[int(item[1])] = item[0]
        
        cDNA_record = FastaRecord( r.name, cDNA_consensus_sequence )
        cDNA_writer.writeRecord( cDNA_record )
        write_fasta( [cDNA_record], temp_fasta )

        ### now that we have established the cDNA sequence ... we start all over again with a cDNA MSA!
        ### get cDNA profile fn
        MSA_cDNA_info_fn = MSA_cDNA_info_fn_dict[locus]
        MSA_cDNA_info = parse_cdna_info( MSA_cDNA_info_fn )

        MSA_cDNA_fn = MSA_cDNA_fn_dict[locus]

        ### add artificial cDNA sequence to the cDNA profile
        output = os.path.join(output_dir, name + "_cDNA.afa")
        add_to_alignment(temp_fasta, MSA_cDNA_fn, output)

        ### read out the same reference sequence from the NEW profile
        ### as well as our artificial cDNA consensus sequence
        reference = get_reference_record( MSA_cDNA_fn )
        aligned_reference = get_aligned_record( output, reference.name )
        aligned_query = get_aligned_record( output, r.name )

        ### create annotation gff3 for this sequence
        ### comparing the same ref between the new and old profile allows us to keep track of the "canonical" coordinates
        #gff3_annotation = create_annotation(comparison_seq, new_comparison_seq, consensus_seq, cDNA_coverage_map, MSA_cDNA_info, r.name, cDNA_quality_map, coord_dict = cDNA_pos_to_gDNA_pos )  
        gff3_annotation = create_annotation2( reference.sequence, 
                                              aligned_reference.sequence, 
                                              aligned_query.sequence, 
                                              MSA_cDNA_info, 
                                              r.name, 
                                              coord_dict = cDNA_pos_to_gDNA_pos )  
        for record in gff3_annotation:
            cDNA_annot_writer.writeRecord(record)
        ## find best match among all profiles for our consensus sequence    
        best_cDNA, match, total, score, cDNA_var_map = MSA_aligner( aligned_query.sequence, 
                                                                    output, 
                                                                    r.name ) 
        best_cDNA = trim_cdna_type( best_cDNA )
        cDNA_allele_handle.write("%s %s %s %s %s %s\n" % ( locus, 
                                                           r.name, 
                                                           best_cDNA, 
                                                           total-match, 
                                                           total, 
                                                           100*float(score) ))
        os.remove(temp_fasta)
    return
    
def trim_cdna_type( hla_type ):
    parts = hla_type.split(':')
    if len(parts) < 4:
        return hla_type
    return ':'.join(parts[:-1])

def parse_genomic_info( info_file ):
    data = {}
    with open(info_file, "r") as handle:
        for line in handle:
            pos1, pos2, region = line.strip().split()
            data[int(pos1)] = info._make( ( pos2, region, 0 ) )
    return data

def parse_cdna_info( info_file ):
    data = {}
    with open(info_file, "r") as handle:
        for line in handle:
            pos1, pos2, region, pos3 = line.strip().split()
            data[int(pos1)] = info._make( ( pos2, region, pos3 ) )
    return data

@memoize
def get_reference_record( alignment_file ):
    max_length = 0
    for record in FastaReader( alignment_file ):
        length = len(re.findall('[AGCT]', record.sequence))
        if length > max_length:
            max_length = length
            longest_record = record
    return longest_record

def get_aligned_record( alignment_file, name ):
    for record in FastaReader( alignment_file ):
        if record.name == name:
            return record

def add_to_alignment(fasta_file, alignment_file, output_file):
    if os.path.isfile( output_file ):
        log.info('Existing MSA alignment found for "%s", skipping...' % fasta_file)
        return
    log.info('Running MUSCLE to align "%s" to the canonical MSA' % fasta_file)
    muscle_args = { 'profile': True,
                    'in1': alignment_file,
                    'in2': fasta_file,
                    'out': output_file }
    run_muscle( muscle_args )

if __name__ == '__main__':
    import argparse

    logging.basicConfig( level=logging.INFO )

    desc = "A tool for annotating genomic HLA sequences"
    parser = argparse.ArgumentParser( description=desc )

    add = parser.add_argument
    add('msa_fofn',
        metavar='MSA',
        help='A FOFN of Multiple Sequence Alignments and information')
    add('input_file',
        metavar='FASTA',
        help='A FASTA file of sequence data to annotate')
    add('locus_dict',
        metavar='LOCUS',
        help='A file specifying which Locus to use to align each sequence')
    add('-o', '--output',
        metavar='DIR',
        default='annotation',
        help='Directory to output results to')

    args = parser.parse_args()
    type_hla( args.msa_fofn,
              args.input_file,
              args.locus_dict,
              args.output )
